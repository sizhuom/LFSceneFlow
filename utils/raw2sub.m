function [ lf ] = raw2sub( raw, sz )
%RAW2SUB Convert a light field from "raw" format to subaperture format

lf = LF2Sub(raw2LF(raw, sz));

end

