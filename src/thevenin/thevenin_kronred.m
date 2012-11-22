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

function [v, timetotal, timelu, timekronred, timebacksolve] = thevenin_kronred(M,N,Nnvc,Nvc)
% 
% Thevenin impedance calculation with kronred and KLU factorization
%
% M must be ordered non-voltage controlled first
%

cutDegree = 3;

index = 1:N; ops1 = []; ops2 = [];
[Manalyze,Nnvckranalyze,index,ops1,ops2,Mkranalyze] = kronred_analyze(M,Nnvc,cutDegree,[],index,false,ops1,ops2,false);
[Manalyze,Nnvckranalyze,index,ops1,ops2,Mkranalyze] = kronred_analyze(Mkranalyze,Nnvckranalyze,cutDegree,Manalyze,index,false,ops1,ops2,false);
[Manalyze,Nnvckranalyze,index,ops1,ops2,Mkranalyze] = kronred_analyze(Mkranalyze,Nnvckranalyze,cutDegree,Manalyze,index,false,ops1,ops2,false);
filter = zeros(1,N); filter(index) = 1:(Nnvckranalyze+Nvc);
Mkrslow = kronred_slow(Manalyze,ops1,index,filter,true);
[Mkrfast timekronred] = kronred(Manalyze,ops1,index,filter);
% assert(norm(Mkrslow-Mkr,1) < 1e-10);
% assert(norm(Mkrfast-Mkr,1) < 1e-10);
assert(norm(Mkrfast-Mkrslow,1) < 1e-10);

% % stats
% [Nnvc-Nnvckranalyze size(ops1,2) size(ops2,2)]

% KLU part
M = Mkrfast;
N = size(M,2);
Nnvc = Nnvckranalyze;

nvc = 1:Nnvc;
vc = (Nnvc+1):N;
Mnvc = M(nvc,nvc);

if Nnvc == 0
    v = 1./diag(M(vc,vc));
    timelu = 0;
    timebacksolve = 0;
    timetotal = timelu + timekronred + timebacksolve;        
    return;
end

[resklu timelu] = thevenin_klu(Mnvc);
L = resklu.L; U = resklu.U; p = resklu.p; q = resklu.q;
R = resklu.R; F = resklu.F;
assert(norm(F,1) == 0);

Mvc = M(vc,vc);
reachL = thevenin_reach(int64(Nnvc),int64(Nvc),L,M(1:Nnvc,(Nnvc+1):end),p);
reachU = thevenin_reach(int64(Nnvc),int64(Nvc),U',M((Nnvc+1):end,1:Nnvc)',q);
[S vfast timebacksolve timeips] = thevenin_getz(int64(Nnvc),int64(Nvc),L,U',R,M(1:Nnvc,(Nnvc+1):end),M((Nnvc+1):end,1:Nnvc)',Mvc,p,q,reachL,reachU,true,'');
timebacksolve = timebacksolve + timeips;

v = zeros(Nvc,1);

for k = 1:Nvc
    ystark = M(1:Nnvc,Nnvc+k);
    ystark = R\ystark(p);
    ustark = cs_lsolve(L,full(ystark));
    assert(norm(L*ustark-ystark,1) < 1e-10);
    
    ykstar = M(Nnvc+k,1:Nnvc)';
    ykstar = ykstar(q);
    lkstar = cs_utsolve(U,full(ykstar));
    
    uk = M(Nnvc+k,Nnvc+k)-lkstar'*ustark;
    
    v(k) = 1/uk;
end

assert(norm(v-vfast,1) < 1e-10);

timetotal = timelu + timekronred + timebacksolve;
fprintf('thevenin_kronred timing: total %e, kronred %e, lu %e, backsolve %e\n',timetotal,timelu,timekronred,timebacksolve);

figure(12)
spy(L)
hold on
spy(U,'r')
title('Y_{nc} Factorization with Elimination')
xlabel(sprintf('nz = %d',nnz(L+U)))

end
