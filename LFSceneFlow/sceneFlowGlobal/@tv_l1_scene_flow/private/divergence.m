function [ div ] = divergence( p )
%DIVERGENCE Compute the divergence of the dual variable using backward
%difference

padx = padarray(p(:,:,:,:,1), [0 1 0 0], 0, 'pre');
padx(:,end,:,:) = 0;
dx = diff(padx, 1, 2);

pady = padarray(p(:,:,:,:,2), [1 0 0 0], 0, 'pre');
pady(end,:,:,:) = 0;
dy = diff(pady, 1, 1);

padu = padarray(p(:,:,:,:,3), [0 0 0 1], 0, 'pre');
padu(:,:,:,end) = 0;
du = diff(padu, 1, 4);

padv = padarray(p(:,:,:,:,4), [0 0 1 0], 0, 'pre');
padv(:,:,end,:) = 0;
dv = diff(padv, 1, 3);

div = dx + dy + du + dv;

end

