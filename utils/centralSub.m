function [ sa ] = centralSub( lf )
%CENTRALSUB Get the central sub-aperture image of an LF

sz = size(lf);
sa = squeeze(lf(round((1+sz(1))/2),round((1+sz(2))/2),:,:,:));

end

