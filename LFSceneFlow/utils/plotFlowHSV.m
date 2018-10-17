function f = plotFlowHSV( flowx, flowy, flowz )
%PLOTFLOWHSV Plot xy flow in HSV, z flow separately
% H: angle from x axis, counter-clockwise, 0-360
% S: magnitude, normalized by maximum flow
% V: occlusion, not implemented yet

% HSV representation of the x-y components
angle = atan2(flowy, flowx);
angle(angle<0) = angle(angle<0) + 2*pi;
h = angle / (2*pi);

mag = sqrt(flowx.^2+flowy.^2);
mag(isinf(mag)) = NaN;
mag(isnan(mag)) = max(mag(:));
flowz(isinf(flowz)) = NaN;
flowz(isnan(flowz)) = max(flowz(:));
maxMag = max(max(mag(:)),max(abs(flowz(:))));
s = mag / maxMag;

hsv = cat(3,h,s,double(~isnan(h)));
flowrgb = hsv2rgb(hsv);

% HSV color wheel for x/y flow
resolution = 1000;
step = maxMag * 2 / resolution;
[x, y] = meshgrid(-maxMag:step:maxMag);
[theta, rho] = cart2pol(x,y);
theta(theta<0) = theta(theta<0) + 2*pi;
h = theta / (2*pi);
s = min(1, rho / maxMag);
v = rho < maxMag;
wheelrgb = hsv2rgb(cat(3,h,s,v));

% HSV colormap for z flow
cmap = ones(64,1,3);
cmap(1:32,1,1) = 0.5;
cmap(1:32,1,2) = 1:-1/32:1/32;
cmap(33:64,1,1) = 0;
cmap(33:64,1,2) = 0:1/32:31/32;
cmap = hsv2rgb(cmap);
cmap = permute(cmap, [1 3 2]);

% Draw figure
axisColor = 'black';
f = figure; 
subplot(2, 2, 1);
imshow(flowrgb);
set(gca,'XColor',axisColor);
set(gca,'YColor',axisColor);
title('XY Flow','Color',axisColor);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'TickLength',[0 0]);
subplot(2, 2, 2);
iptsetpref('ImshowAxesVisible','on');
imshow(wheelrgb,'XData',[-maxMag,maxMag],'YData',[-maxMag,maxMag]);
ticks = -maxMag:maxMag/2:maxMag;
tickLabels = arrayfun(@(x) sprintf('%.1f',x), ticks, 'UniformOutput', false);
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);
set(gca,'YTick',ticks);
set(gca,'YTickLabel',tickLabels);
set(gca,'XColor',axisColor);
set(gca,'YColor',axisColor);
title('XY Color Coding','Color',axisColor);
subplot(2, 2, [3 4]);
imagesc(flowz,[-maxMag,maxMag]); 
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'TickLength',[0 0]);
colorbar('XColor',axisColor,'YColor',axisColor,'Ticks',ticks,'TickLabels',tickLabels);
% colormap parula;
colormap(cmap);
pbaspect([size(flowx,2) size(flowy,1) 1]);
set(gca,'XColor',axisColor);
set(gca,'YColor',axisColor);
title('Z Flow','Color',axisColor);

end

