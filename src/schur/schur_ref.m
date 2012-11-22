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

function [S, timetotal, timelu, timekronred, timebacksolve] = thevenin_ref(M,N,Nnvc,Nvc)
% 
% reference Schur complement calculation
%
% M must be ordered non-voltage controlled first
%

nvc = 1:Nnvc;
vc = (Nnvc+1):N;
Mnvc = M(nvc,nvc);
Mvc = M(vc,vc);

S = Mvc-M((Nnvc+1):end,1:Nnvc)*(Mnvc\M(1:Nnvc,(Nnvc+1):end));

timetotal = 0.0; % no timing here
timelu = 0.0;
timekronred = 0.0;
timebacksolve = 0.0;


end
