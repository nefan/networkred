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

function [v, timetotal, timelu, timekronred, timebacksolve] = thevenin_ref(M,N,Nnvc,Nvc)
% 
% factor-solve Thevenin impedance calculation
%
% M must be ordered non-voltage controlled first
%

nvc = 1:Nnvc;
vc = (Nnvc+1):N;
Mnvc = M(nvc,nvc);

[L,U,p,q] = cs_lu(Mnvc,0.0);

v = zeros(Nvc,1);

for k = 1:Nvc
    ystark = M(1:Nnvc,Nnvc+k);
    ystark = ystark(p);
    ustark = cs_lsolve(L,full(ystark));
    assert(norm(L*ustark-ystark,1) < 1e-10);
    
    ykstar = M(Nnvc+k,1:Nnvc)';
    ykstar = ykstar(q);
    lkstar = cs_utsolve(U,full(ykstar));
    
    uk = M(Nnvc+k,Nnvc+k)-lkstar'*ustark;
    
    v(k) = 1/uk;
end

timetotal = 0.0; % no timing here
timelu = 0.0;
timekronred = 0.0;
timebacksolve = 0.0;


end
