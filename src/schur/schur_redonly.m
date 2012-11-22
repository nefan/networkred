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
% Schur complement calculation with node eliminaton
%
% M must be ordered non-voltage controlled first
%

cutDegree = N+1;

nvc = 1:Nnvc;
vc = (Nnvc+1):N;
Mnvc = M(nvc,nvc);
Mvc = M(vc,vc);

Pnvc = amd(Mnvc);
Pvc = 1:Nvc;
Phj = [Pnvc (Nnvc+Pvc)];
M = M(Phj,Phj);

index = 1:N; ops1 = []; ops2 = [];
[Manalyze,Nnvckranalyze,index,ops1,ops2,Mkranalyze] = kronred_analyze(M,Nnvc,cutDegree,[],index,false,ops1,ops2,true);
filter = zeros(1,N); filter(index) = 1:(Nnvckranalyze+Nvc);
% Mkrslow = kronred_slow(Manalyze,ops1,index,filter,true);
[Mkrfast timekronred] = kronred(Manalyze,ops1,index,filter);
% assert(norm(Mkrfast-Mkrslow,1) < 1e-10);

% KLU part
M = Mkrfast;
N = size(M,2);
Nnvc = Nnvckranalyze;

assert(Nnvc == 0);
S = M;
timelu = 0;
timebacksolve = 0;
timeips = 0;
timetotal = timelu + timekronred + timebacksolve + timeips;        

fprintf('schur_redonly timing: total %e, kronred %e\n',timetotal,timekronred);

figure(13)
spy(M)
title('M reduced')


end
