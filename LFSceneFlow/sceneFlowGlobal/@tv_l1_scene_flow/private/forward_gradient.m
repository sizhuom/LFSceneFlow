function [ gx, gy, gu, gv ] = forward_gradient( v )
%FORWARD_GRADIENT Compute the forward gradient of v

gx = padarray(diff(v,1,2), [0 1 0 0], 0, 'post');
gy = padarray(diff(v,1,1), [1 0 0 0], 0, 'post');
gu = padarray(diff(v,1,4), [0 0 0 1], 0, 'post');
gv = padarray(diff(v,1,3), [0 0 1 0], 0, 'post');

end

