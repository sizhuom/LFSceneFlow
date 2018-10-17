function [ trans ] = buildMotion( qw )
%BUILDMOTION Build the transformation matrix from the qw vector
qw = reshape(qw, [6 1]);
trans = [rodrigues(qw(4:6)) qw(1:3);
         zeros(1, 3) 1];

end

