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

function [S, timetotal, timelu, timekronred, timebacksolve] = schur_klu(M,N,Nnvc,Nvc,name)
% 
% Schur complement calculation with KLU factorization and kron reduction
%
% M must be ordered non-voltage controlled first
%

cutDegree = 3;

index = 1:N; ops1 = []; ops2 = [];
[Manalyze,Nnvckranalyze,index,ops1,ops2,Mkranalyze] = kronred_analyze(M,Nnvc,cutDegree,[],index,false,ops1,ops2,true);
[Manalyze,Nnvckranalyze,index,ops1,ops2,Mkranalyze] = kronred_analyze(Mkranalyze,Nnvckranalyze,cutDegree,Manalyze,index,false,ops1,ops2,true);
[Manalyze,Nnvckranalyze,index,ops1,ops2,Mkranalyze] = kronred_analyze(Mkranalyze,Nnvckranalyze,cutDegree,Manalyze,index,false,ops1,ops2,true);
filter = zeros(1,N); filter(index) = 1:(Nnvckranalyze+Nvc);
Mkrslow = kronred_slow(Manalyze,ops1,index,filter,true);
[Mkrfast timekronred] = kronred(Manalyze,ops1,index,filter);
assert(norm(Mkrfast-Mkrslow,1) < 1e-10);

% KLU part
M = Mkrfast;
N = size(M,2);
Nnvc = Nnvckranalyze;

nvc = 1:Nnvc;
vc = (Nnvc+1):N;
Mnvc = M(nvc,nvc);
Mvc = M(vc,vc);

if Nnvc == 0
    S = Mvc;
    timelu = 0;
    timebacksolve = 0;
    timeips = 0;
    timetotal = timelu + timekronred + timebacksolve + timeips;        
    return;
end

[resklu timelu] = thevenin_klu(Mnvc);
L = resklu.L; U = resklu.U; p = resklu.p; q = resklu.q;
R = resklu.R; F = resklu.F;
assert(norm(F,1) == 0);

reachL = thevenin_reach(int64(Nnvc),int64(Nvc),L,M(1:Nnvc,(Nnvc+1):end),p);
reachU = thevenin_reach(int64(Nnvc),int64(Nvc),U',M((Nnvc+1):end,1:Nnvc)',q);

[S vfast timebacksolve timeips] = thevenin_getz(int64(Nnvc),int64(Nvc),L,U',R,M(1:Nnvc,(Nnvc+1):end),M((Nnvc+1):end,1:Nnvc)',Mvc,p,q,reachL,reachU,false,name);

timetotal = timelu + timekronred + timebacksolve + timeips;
fprintf('schur_kronred timing: total %e, kronred %e, lu %e, backsolve %e, ips %e\n',timetotal,timekronred,timelu,timebacksolve,timeips);
if exist(['../tmp/' name '.from-cuda'],'file')
    fromcuda = importdata(['../tmp/' name '.from-cuda']);
    timecuda = fromcuda(1);
    fprintf('schur_kronred timing: equivalent CUDA ips %e (error %e)\n',timecuda,fromcuda(2));
end

figure(11)
spy(L)
hold on
spy(U,'r')
title('Y_{nc} Factorization')

end
