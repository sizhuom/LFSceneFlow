function [ w ] = makeOmega( degree, axis )
%MAKEOMEGA Make a rotation vector
rad = degree / 180 * pi;
if (ischar(axis))
    switch axis
        case 'x'
            w = [rad; 0; 0];
        case 'y'
            w = [0; rad; 0];
        case 'z'
            w = [0; 0; rad];
    end
else
    w = rad * axis;
end

end
