function [ raw ] = LF2Raw( lf )
%LF2RAW Convert a light field from 4D array format to "raw" format

raw = permute(lf, [1 3 2 4 5]);
raw = reshape(raw, size(lf,1)*size(lf,3), size(lf,2)*size(lf,4), size(lf,5));


end

