%
% This file is part of networkred.
%
% Copyright (C) 2012, Technical University of Denmark
% https://github.com/nefan/networkred
%
% networkred is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% networkred is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with networkred.  If not, see <http://www.gnu.org/licenses/>.
% 

function [Manalyze, Nnvckr, index, ops1, ops2, Mkr] = kronred_analyze(M,Nnvc,cutDegree,Manalyze,oldindex,split,ops1,ops2,updateEntireMatrix)

Mkr = M;
Nnvckr = Nnvc;
N = size(M,2);

index = oldindex;

krI = 1; % index into Mkr
% stat = zeros(3,10000);
stat = [0 0 0];
for i=1:Nnvc    
    node = oldindex(i);
    col = Mkr(:,krI);
    degree = nnz(col(1:Nnvckr))-1; % discarding self-loops
%     degree = nnz(col)-1; % discarding self-loops

    if degree <= cutDegree % cut
        
        filter = zeros(size(col));
        filter(krI) = 1;
        
        % update
        row = Mkr(krI,:);
        update = col*row/col(krI);
%         if (nnz(update) > 120)
%             krI = krI + 1;            
%             continue;
%         end        
        Mkr = Mkr - update;        

        Mkr = Mkr(~filter,~filter); % remove node
        index = [index(1:krI-1) index(krI+1:end)];        
        Nnvckr = Nnvckr - 1;   
        
        [ies,jes,xx] = find(update(~filter,~filter));
        nr1 = 0; nr2 = 0;
        lops1 = [0 node]'; lops2 = [0 node]';
        for ke = 1:length(ies)
            ie = ies(ke); je = jes(ke);
            if (~updateEntireMatrix && min(ie,je) > Nnvckr && ie ~= je)
                continue;
            end            
            if (~split || max(ie,je) <= Nnvckr)
                lops1 = [lops1 [index(ie) index(je)]'];            
                nr1 = nr1 + 1;
            else
                lops2 = [lops2 [index(ie) index(je)]'];            
                nr2 = nr2 + 1;
            end                                    
        end
        while (mod(nr1,4) ~= 0)
            lops1 = [lops1 [0 -1]'];            
            nr1 = nr1 + 1;
        end
        while (mod(nr2,4) ~= 0)
            lops2 = [lops2 [0 -1]'];            
            nr2 = nr2 + 1;
        end        
        assert(mod(nr1,4) == 0); % simd parallelization
        assert(mod(nr2,4) == 0); % simd parallelization
        lops1(1,1) = -nr1;
        lops2(1,1) = -nr2;
        ops1 = [ops1 lops1]; % must be run
%         if nr2 > 0
            ops2 = [ops2 lops2];        
%         end
    else
        krI = krI + 1;
    end
    
end

if isempty(Manalyze)
    Manalyze = M;
end
ops = [ops1 ops2];
for i = 1:size(ops,2)
    op = ops(:,i)';
    if (op(1) <= 0)
        node1 = op(2);
        node2 = op(2);
    else
        node1 = op(1);
        node2 = op(2);
    end    
    if (node2 > 0 && Manalyze(node1,node2) == 0)
        Manalyze(node1,node2) = nan;
    end
end

end
