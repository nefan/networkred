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

function Mkr = kronred_matlab(Msymbolic,ops,index,filter,doRemove)

Mkr = Msymbolic;

i = 1;
while (i <= size(ops,2))
    op = ops(:,i);
    
    assert(op(1) <= 0);
    
    K = -op(1);
    node = op(2);
    nodev = Mkr(node,node);
%     assert(mod(K,4) == 0); % simd
    i = i + 1;
    for j=1:K
        op = ops(:,i);
        if (op(2) > 0)
            assert(Mkr(op(1),op(2)) ~= 0);
            if (isnan(Mkr(op(1),op(2))))
                Mkr(op(1),op(2)) = 0;
            end
            Mkr(op(1),op(2)) = Mkr(op(1),op(2)) - Mkr(op(1),node)*Mkr(node,op(2))/nodev;
        end

        i = i + 1;
    end    
end

if doRemove
    Mkr = Mkr(index,index); % remove discarded nodes
end

end
