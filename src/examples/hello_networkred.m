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

% Example Kron reduction
% This MATLAB file mimics the code in hello_networkred.cpp

addpath('../../lib');

M = sparse(complex([1 .1 0 0.4; .1 2 .2 0; 0 .2 3 .3; .4 0 .3 4]));
Nnvc = 1; % one node to be removed

cutDegree = 4; % reduce all nvc nodes
Nsweeps = 1; % only one sweep necessary with this cutDegree

NITER = 1; % only for performance testing

[Msymbolic, Nkr, index, ops1, ops2] = kronred_symbolic(M,int64(Nnvc),int64(cutDegree),int64(Nsweeps),false,true);
[Mkr timekronred] = kronred(Msymbolic,ops1,index,int64(NITER));

fprintf('M before reduction:\n');
full(M)
fprintf('M after reduction:\n');
full(Mkr)