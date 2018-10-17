function f = plotFlowHSV2D( flowx, flowy )
%PLOTFLOWHSV2D Plot xy flow in HSV
% H: angle from x axis, counter-clockwise, 0-360
% S: magnitude, normalized by maximum flow
% V: occlusion, not implemented yet

% HSV representation of the x-y components
angle = atan2(flowy, flowx);
angle(angle<0) = angle(angle<0) + 2*pi;
h = angle / (2*pi);

mag = sqrt(flowx.^2+flowy.^2);
maxMag = max(mag(:));
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

% Draw figure
axisColor = 'black';
f = figure; 
subplot(1, 2, 1);
imshow(flowrgb);
set(gca,'XColor',axisColor);
set(gca,'YColor',axisColor);
title('XY Flow','Color',axisColor);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'TickLength',[0 0]);
subplot(1, 2, 2);
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

end

