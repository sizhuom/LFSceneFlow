%% compute alpha
global LFTopDir
workDir = fullfile(LFTopDir, 'Images/Sim/04-11-CVPR18', '0714-sceneflow');
ccFilename = 'cc-frame-480_0000.png';
frameNum = 0;
gtAlpha = calcGTAlpha(workDir, ccFilename, frameNum);

%% visualize central subaperture
gtAlphaC = squeeze(gtAlpha(8,8,:,:));
figure; imagesc(gtAlphaC); colorbar;

%% interactive EPI visualization
lf = imread(fullfile(workDir, 'frame-480_0000.png'));
lf = raw2LF(lf, [15 15 240 320]);
y = 8; v = 120;
epi = squeeze(lf(y,:,v,:,:));
u = 160; x = 8;
f1 = figure;
imagesc(epi);
pbaspect([320 15 1]);
hold on;
plot([u-gtAlpha(y,x,v,u)*(1-x),u-gtAlpha(y,x,v,u)*(15-x)],[1 15], '-r');
hold off;
gtAlpha(y,x,v,u)
[u-gtAlpha(y,x,v,u)*(1-x),u-gtAlpha(y,x,v,u)*(15-x)]
while true
    figure(f1);
    [u,x] = ginput(1);
    u = round(u);
    if u < 1
        u = 1;
    elseif u > 320
        u = 320;
    end
    x = round(x);
    if x < 1
        x = 1;
    elseif x > 15
        x = 15;
    end
    
    gtAlpha(y,x,v,u)
    [u-gtAlpha(y,x,v,u)*(1-x),u-gtAlpha(y,x,v,u)*(15-x)]
    
    imagesc(epi);
    pbaspect([320 15 1]);
    hold on;
    plot([u-gtAlpha(y,x,v,u)*(1-x),u-gtAlpha(y,x,v,u)*(15-x)],[1 15], '-r');
    hold off;
    
end