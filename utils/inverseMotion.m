function [ iqw ] = inverseMotion( qw )
%INVERSEMOTION compute the inverse of the motion defined by [q;w]
T = qw(1:3);
invR = eye(3) + crossMatrix(-qw(4:6));
iqw = [-invR*T;-qw(4:6)];

end

