function [ gx, gy, gu, gv ] = backward_gradient( v )
%BACKWARD_GRADIENT Compute the backward gradient of v

padx = padarray(v, [0 1 0 0], 0, 'pre');
padx(:,end,:,:) = 0;
gx = diff(padx, 1, 2);

pady = padarray(v, [1 0 0 0], 0, 'pre');
pady(end,:,:,:) = 0;
gy = diff(pady, 1, 1);

padu = padarray(v, [0 0 0 1], 0, 'pre');
padu(:,:,:,end) = 0;
gu = diff(padu, 1, 4);

padv = padarray(v, [0 0 1 0], 0, 'pre');
padv(:,:,end,:) = 0;
gv = diff(padv, 1, 3);

end

