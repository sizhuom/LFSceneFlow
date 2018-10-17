function [ R ] = rodrigues( w )
%RODRIGUES Make a rotation matrix using the Rodrigues' formula

theta = norm(w);
if theta == 0
    R = eye(3);
else
    w = w / theta;
    wx = crossMatrix(w);
    if theta < 1e-2
        R = eye(3) + theta * wx;
    else
        R = eye(3) + sin(theta)*wx + (1-cos(theta))*wx^2;
    end
end

end

