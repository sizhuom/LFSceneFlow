function [ R, t ] = buildLookAt( pos, lookAt, up, backward )
%BUILDLOOKAT build a lookAt matrix
% forward: world -> camera
% backward: camera -> world

if nargin < 3
    up = [0; 1; 0];
end

z = lookAt(:) - pos(:);
z = z / norm(z);
x = cross(up(:), z);
x = x / norm(x);
y = cross(z, x);

if nargin >= 4 && backward
    R = [x y z];
    t = pos(:);
else
    R = [x y z]';
    t = -R * pos(:);
end

end

