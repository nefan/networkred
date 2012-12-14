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

% kron reduction without symbolic step. For testing
function [Mkr,Nnvckr] = kreduce(M,Nnvc,cutDegree)

Mkr = M;
Nnvckr = Nnvc;

krI = 1; % index into Mkr
for i=1:Nnvc    
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
        Nnvckr = Nnvckr - 1;        
    else
        krI = krI + 1;
    end
    
end


end
