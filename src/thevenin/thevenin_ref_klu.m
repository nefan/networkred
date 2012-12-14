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

function [v, timetotal, timelu, timekronred, timebacksolve] = thevenin_ref_klu(M,N,Nnvc,Nvc)
% 
% semi smart Thevenin impedance calculation with KLU factorization
%
% M must be ordered non-voltage controlled first
%

global NITER;

nvc = 1:Nnvc;
vc = (Nnvc+1):N;
Mnvc = M(nvc,nvc);
Mvc = M(vc,vc);

[resklu timelu] = thevenin_klu(Mnvc,NITER);
L = resklu.L; U = resklu.U; p = resklu.p; q = resklu.q;
R = resklu.R; F = resklu.F;
assert(norm(F,1) == 0);

reachL = thevenin_reach(int64(Nnvc),int64(Nvc),L,M(1:Nnvc,(Nnvc+1):end),p);
reachU = thevenin_reach(int64(Nnvc),int64(Nvc),U',M((Nnvc+1):end,1:Nnvc)',q);

timekronred = 0;

[S vfast  timebacksolve timeips] = thevenin_getz(int64(Nnvc),int64(Nvc),L,U',R,M(1:Nnvc,(Nnvc+1):end),M((Nnvc+1):end,1:Nnvc)',Mvc,p,q,reachL,reachU,true,'',NITER);
timebacksolve = timebacksolve + timeips;
% norm(S-(Mvc-M((Nnvc+1):end,1:Nnvc)*(Mnvc\M(1:Nnvc,(Nnvc+1):end))),1)

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
fprintf('thevenin_ref_klu timing: total %e, lu %e, backsolve %e\n',timetotal,timelu,timebacksolve);

figure(11)
spy(L)
hold on
spy(U,'r')
title('Y_{nc} Factorization')

end
