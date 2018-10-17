function [ sub ] = LF2Sub( lf )
%LF2SUB Convert a light field from 4D array format to subaperture format

sub = permute(lf, [3 1 4 2 5]);
sub = reshape(sub, size(lf,1)*size(lf,3), size(lf,2)*size(lf,4), size(lf,5));

end

