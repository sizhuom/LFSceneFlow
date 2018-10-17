function [ lf ] = raw2LF( raw, sz )
%RAW2LF Convert a light field from "raw" format to 4D array format

lf = reshape(raw, sz(1), sz(3), sz(2), sz(4), size(raw, 3));
lf = permute(lf, [1 3 2 4 5]);


end

