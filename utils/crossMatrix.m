function [ ax ] = crossMatrix( a )
%CROSSMATRIX Cross product matrix of a 3-vector
%   crossMatrix(a)*b == cross(a,b)
ax = [0 -a(3) a(2);
      a(3) 0 -a(1);
      -a(2) a(1) 0];

end

