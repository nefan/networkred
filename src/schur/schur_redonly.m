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

global NITER;

cutDegree = N+1;

nvc = 1:Nnvc;
vc = (Nnvc+1):N;

[Msymbolic, Nkr, index, ops1, ops2] = kronred_symbolic(M,int64(Nnvc),int64(cutDegree),int64(3),false,true);
[Mkr timekronred] = kronred(Msymbolic,ops1,index,NITER);

% KLU part
M = Mkr;
N = size(M,2);
Nnvc = Nkr-Nvc;

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
