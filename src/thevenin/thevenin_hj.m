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

function [v, timetotal, timelu, timekronred, timebacksolve] = thevenin_hj(M,N,Nnvc,Nvc)
% 
% Thevenin impedance calculation by Johansson et al.
%
% M must be ordered non-voltage controlled first
%

nvc = 1:Nnvc;
vc = (Nnvc+1):N;
Mnvc = M(nvc,nvc);
Mvc = M(vc,vc);

Pnvc = amd(Mnvc);
Pvc = amd(Mvc);
Phj = [Pnvc (Nnvc+Pvc)];
[L,U,p,q,timelu] = thevenin_umfpack(M(Phj,Phj),0);

v = zeros(Nvc,1);

for k = 1:Nvc
    ustark = U(1:Nnvc,Nnvc+k);    
    lkstar = L(Nnvc+k,1:Nnvc);
    
    uk = M(Phj(Nnvc+k),Phj(Nnvc+k))-lkstar*ustark;
    
    v(Pvc(k)) = 1/uk;
end

timekronred = 0;
timebacksolve = 0;
timetotal = timelu + timekronred + timebacksolve;
fprintf('thevenin_umfpack timing: total %e, lu %e, backsolve %e\n',timetotal,timelu,timebacksolve);

figure(10)
spy(L)
hold on
spy(U,'r')
hv=line([Nnvc,Nnvc],[0 N]);
hh=line([0 N],[Nnvc,Nnvc]);
set(hv,'LineStyle','--','LineWidth',3,'Color','k');
set(hh,'LineStyle','--','LineWidth',3,'Color','k');
title('Y Factorization')
xlabel(sprintf('nz = %d',nnz(L+U)))

end
