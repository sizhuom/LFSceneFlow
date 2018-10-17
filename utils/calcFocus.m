function [ fLen ] = calcFocus( camParam, focus )
%CALCFOCUS Calculate the correct focal length to focus at a give depth

[~, param] = genIntrinsicsP2(camParam.resol, camParam.apert, camParam.fov, 1, [0 0]);
fLen = 1/(1/(param.dm+param.dM)+1/focus);

end

