/*
 * This file is part of networkred.
 *
 * Copyright (C) 2012, Technical University of Denmark
 * https://github.com/nefan/networkred
 *
 * networkred is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * networkred is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with networkred.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "kronred_symbolic.hpp"
//#undef NDEBUG
//#define assert(X) mxAssert(X,"c assert")

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
using namespace boost;

// graph
typedef property<vertex_index_t, int,
        property<vertex_color_t, bool> > Vertex_t;
typedef property<edge_index_t, int,
        property<edge_weight_t,std::complex<double> > > Edge_t;
typedef subgraph<adjacency_list<setS, vecS, bidirectionalS, Vertex_t, Edge_t> > Graph;

typedef typename graph_traits<Graph>::edge_descriptor Edge;
typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::vertex_iterator vertex_iter;
typedef graph_traits<Graph>::in_edge_iterator in_edge_iter;
typedef graph_traits<Graph>::out_edge_iterator out_edge_iter;

typedef property_map<Graph, vertex_index_t>::type vertex_indicies;
typedef property_map<Graph, vertex_color_t>::type colors;
typedef property_map<Graph, edge_index_t>::type edge_indicies;
typedef property_map<Graph, edge_weight_t>::type weights;

// number of attached active nodes
int active_out_degree(Vertex v, Graph& g) { 
    int out_degree = 0;
    // loop through pairs of edges
    std::pair<out_edge_iter, out_edge_iter> eop;
    for (eop = out_edges(v,g); eop.first != eop.second; ++eop.first) {
        Vertex w = target(*eop.first, g);
        if (get(vertex_color, g, w)) // active node
            out_degree++;
    }

    return out_degree;
}

// encode indices corresponding to M into operations
void kronred_encode(const cljmexComplex_sparse_matrix& M, const Ops& ops, EncodedOps& eOps)
{
    for (Ops::const_iterator oit = ops.begin(); oit != ops.end();) {
        const Op& op = *oit;

        int K = -op.first;
        assert(K >= 0); // ops node indexing failure
#ifdef USE_AVX_SIMD
        assert(K % 4 == 0); // ops node indexing failure (fail mod 4)
#endif

        int node = op.second;
        assert(node >= 0); // ops indexing failure

        {
            EncodedOp eop;
            eop[0] = -K;
            eop[1] = linear_index(node,node,M.p,M.i);
            eop[2] = 0;
            eop[3] = 0;

            eOps.push_back(eop);
        }

        oit++;
        for (int j=0; j<K; j +=1, oit++) {
            EncodedOp eop;

            if (oit->second < 0) { // trash
                eop[0] = M.nz;
                eop[1] = M.nz;
                eop[2] = M.nz;
                eop[3] = 0;
            } else {
                eop[0] = linear_index(oit->first,oit->second,M.p,M.i);
                eop[1] = linear_index(oit->first,node,M.p,M.i);
                eop[2] = linear_index(node,oit->second,M.p,M.i);
                eop[3] = 0;
            }

            eOps.push_back(eop);
        }
    }

    assert(ops.size() == eOps.size());
}

void kronred_symbolic(const cljmexComplex_sparse_matrix& M, 
        const int Nnvc, const int cutDegree, const int Nsweeps, 
        const bool split, const bool updateEntireMatrix,
        cljmexComplex_sparse_matrix& Msymbolic,
        int& Nkr, std::vector<int>& index,
        EncodedOps& eOps1, EncodedOps& eOps2) {

    const int n = M.rows;

    // create graph representing admittance matrix M
    Graph g(n);;
    Graph gnvc = g.create_subgraph();

    // add nodes
    {
        std::pair<vertex_iter, vertex_iter> vp;
        int i;
        for (vp = vertices(g), i = 0; vp.first != vp.second; ++vp.first, i++) {
            Vertex v = *vp.first;
            assert(get(vertex_index, g, v) == i);
            put(vertex_color, g, v, true);
            if (i < Nnvc) // nvc node
                add_vertex(v,gnvc);
        }
    }

    // add edges
    for (int c=0; c<M.cols; c++) // columns
        for (int ri=M.p[c]; ri<M.p[c+1]; ri++) {
            int r = M.i[ri]; // row

            Edge e = add_edge(c,r,g).first;
            put(edge_weight, g, e, M.v[ri]);
            if (std::max(c,r) < Nnvc)
                add_edge(c,r,gnvc);
        }

    // compute operations
    Ops ops1, ops2;
    Nkr = n;
    for (int sweep=0; sweep<Nsweeps; sweep++) { // sweep several times over the nvc nodes
        std::pair<vertex_iter, vertex_iter> vp = vertices(gnvc);
        vertex_iter next;
        for (next = vp.first; vp.first != vp.second; vp.first = next) { // loop over nvc nodes
            ++next;
            Vertex v = *vp.first;
            const int index_node = get(vertex_index, g, v);
            if (!get(vertex_color, g, v)) // inactive node
                continue;

            // active (not already removed) and of low degree
            if (active_out_degree(v,gnvc)-1 <= cutDegree) { // discarding self-loops
                put(vertex_color, g, v, false); // mark inactive
                Nkr--;

                // start operations
                ops1.push_back(StartNodeOp(index_node));
                ops2.push_back(StartNodeOp(index_node));
                int startOp1 = ops1.size()-1;
                int startOp2 = ops2.size()-1;
                int nodeOps1 = 0;
                int nodeOps2 = 0;
                // loop through pairs of edges
                std::pair<out_edge_iter, out_edge_iter> eop;
                for (eop = out_edges(v,g); eop.first != eop.second; ++eop.first) {
                    Vertex wi = target(*eop.first, g);
                    if (!get(vertex_color, g, wi)) // inactive node
                        continue; 

                    const int index_i = get(vertex_index, g, wi);

                    std::pair<in_edge_iter, in_edge_iter> eip;
                    for (eip = in_edges(v,g); eip.first != eip.second; ++eip.first) {
                        Vertex wj = source(*eip.first, g);
                        if (!get(vertex_color, g, wj)) // inactive node
                            continue; 

                        const int index_j = get(vertex_index, g, wj);

                        if (!updateEntireMatrix && std::min(index_i,index_j) >= Nnvc && index_i != index_j) // skip
                            continue;

                        if (!edge(wi,wj,g).second) { // if no already existing edge
                            Edge e = add_edge(wi,wj,g).first; // add new edge to graph
                            put(edge_weight, g, e, std::complex<double>(0,0));
                        }
                        if (!split || std::max(index_i,index_j) < Nnvc) {
                            // internal node (if split) to be added to ops1
                            if (std::max(index_i,index_j) < Nnvc)
                                add_edge(wi,wj,gnvc); // add new edge to nvc subgraph
                            ops1.push_back(NewOp(index_i,index_j));
                            nodeOps1++;
                        } else {
                            // external node (if split) to be added to ops2
                            ops2.push_back(NewOp(index_i,index_j));
                            nodeOps2++;
                        }

                    }
                }

#ifdef USE_AVX_SIMD
                // insert NoOp to make sure ops* are multiple of 4
                while (nodeOps1 % 4 != 0) {
                    ops1.push_back(NoOp);
                    nodeOps1++;
                }
                while (nodeOps2 % 4 != 0) {
                    ops2.push_back(NoOp);
                    nodeOps2++;
                }
#endif
                ops1[startOp1].first = -nodeOps1;
                ops2[startOp2].first = -nodeOps2;
            }
        }
    }

    // output
    Msymbolic.rows = n;
    Msymbolic.cols = n;
    Msymbolic.nz = num_edges(g);
    Msymbolic.p = new cljmexInt[n+1];
    Msymbolic.i = new cljmexInt[Msymbolic.nz];
    Msymbolic.x = new cljmexEntry[Msymbolic.nz];
    Msymbolic.z = new cljmexEntry[Msymbolic.nz];
    {
        std::pair<vertex_iter, vertex_iter> vp = vertices(g);
        int ri = 0;
        for (int c=0; vp.first != vp.second; vp.first++, c++) { // loop over vertices/columns
            Vertex v = *vp.first;
            Msymbolic.p[c] = ri;

            std::pair<out_edge_iter, out_edge_iter> eop;
            for (eop = out_edges(v,g); eop.first != eop.second; eop.first++, ri++) {
                Vertex w = target(*eop.first, g);
                const int r = get(vertex_index, g, w);

                Msymbolic.i[ri] = r;
                std::complex<double> weight = get(edge_weight, g, *eop.first);
                Msymbolic.x[ri] = real(weight);
                Msymbolic.z[ri] = imag(weight);
            }
        }
        assert(ri == num_edges(g));
        Msymbolic.p[n] = ri;
    }

    for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; vp.first++) { // loop over vertices/columns
        Vertex v = *vp.first;

        if (get(vertex_color, g, v)) // active node
            index.push_back(get(vertex_index, g, v));
    }

    kronred_encode(Msymbolic,ops1,eOps1);
    kronred_encode(Msymbolic,ops2,eOps2);
}

cljmex_start()

    const int n = M.rows;

    // symbolic output
    cljmexComplex_sparse_matrix Msymbolic;
    int Nkr; 
    std::vector<int> index;
    EncodedOps eOps1, eOps2;
    Ops ops1;

    // perform symbolic analysis
    kronred_symbolic(M, Nnvc, cutDegree, Nsweeps, 
        split, updateEntireMatrix,
        Msymbolic, Nkr, index, eOps1, eOps2);

    // output
    Msymbolic_out = mxCreateSparse(Msymbolic.rows,Msymbolic.cols,Msymbolic.nz, mxCOMPLEX) ;
    cljmexInt* Msymbolic_p = (cljmexInt *) mxGetJc (Msymbolic_out);
    cljmexInt* Msymbolic_i = (cljmexInt *) mxGetIr (Msymbolic_out);
    double* Msymbolic_x = mxGetPr (Msymbolic_out);
    double* Msymbolic_z = mxGetPi (Msymbolic_out);
    for (int c=0; c<Msymbolic.cols+1; c++)
        Msymbolic_p[c] = Msymbolic.p[c];
    for (int i=0; i<Msymbolic.nz; i++) {
        Msymbolic_i[i] = Msymbolic.i[i];
        Msymbolic_x[i] = Msymbolic.x[i];
        Msymbolic_z[i] = Msymbolic.z[i];
    }

    index_out = mxCreateDoubleMatrix(1,Nkr,mxREAL);
    {
        double* p = (double*)mxGetPr(index_out);
        std::vector<int>::iterator it = index.begin();
        while (it != index.end()) *p++ = (double)(*it++) +1; // +1 for matlab indexing
     }
 
    ops1_out = mxCreateDoubleMatrix(4,eOps1.size(),mxREAL);
    {
        double* p = (double*)mxGetPr(ops1_out);
        EncodedOps::iterator it = eOps1.begin();
        while (it != eOps1.end()) {
            EncodedOp& eop = *it++;
            *p++ = eop[0]; *p++ = eop[1];
            *p++ = eop[2]; *p++ = eop[3];
        }
    }
 
    ops2_out = mxCreateDoubleMatrix(4,eOps2.size(),mxREAL);
    {
        double* p = (double*)mxGetPr(ops2_out);
        EncodedOps::iterator it = eOps2.begin();
        while (it != eOps2.end()) {
            EncodedOp& eop = *it++;
            *p++ = eop[0]; *p++ = eop[1];
            *p++ = eop[2]; *p++ = eop[3];
        }
    }

    // clean up
    delete Msymbolic.p;
    delete Msymbolic.i;
    delete Msymbolic.x;
    delete Msymbolic.z;
 
cljmex_end()

