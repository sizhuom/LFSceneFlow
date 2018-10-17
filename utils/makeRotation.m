function [ R ] = makeRotation( degree, axis )
%MAKEROTATION Make a rotation matrix
rad = degree / 180 * pi;
if (ischar(axis))
    switch axis
        case 'x'
            R = [1 0 0; 0 cos(rad) -sin(rad); 0 sin(rad) cos(rad)];
        case 'y'
            R = [cos(rad) 0 sin(rad); 0 1 0; -sin(rad) 0 cos(rad)];
        case 'z'
            R = [cos(rad) -sin(rad) 0; sin(rad) cos(rad) 0; 0 0 1];
    end
else
    R = rodrigues(rad*axis);
end

end

