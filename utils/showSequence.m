function showSequence( images )
%SHOWSEQUENCE Visualize a sequence of images
curr = 1;
fig = figure('keypressfcn', @(H,E) keypress_callback(H,E,images));
imshow(images{curr});


% Callback
function keypress_callback(~, E, images)
switch E.Key
    case 'rightarrow'
        curr = min(length(images), curr+1);
    case 'leftarrow'
        curr = max(1, curr-1);
end
imshow(images{curr});
end

end