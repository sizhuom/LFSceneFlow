function [ patch, Hnew ] = makePatch( lf, coord, H )
%MAKEPATCH cut a patch out of the original LF

patch = lf(coord(1,1):coord(1,2),coord(2,1):coord(2,2),coord(3,1):coord(3,2),...
    coord(4,1):coord(4,2),:);

if nargin == 3
    offset = eye(5);
    offset(1:4,5) = coord([2 1 4 3],1) - 1;
    Hnew = H * offset;
else
    Hnew = [];
end

end

