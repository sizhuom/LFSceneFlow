function [ L ] = ndlaplacian( n )
%NDLAPLACIAN Make a n-dimensional Laplacian operator

sz = ones(1, n) * 3;
if n == 1
    L = zeros(1, 3);
else
    L = zeros(sz);
end

sub = ones(n, 1) * 2;
L(sub2indLP(sub)) = -2*n;
for i = 1:n
    sub(i) = 1;
    L(sub2indLP(sub)) = 1;
    sub(i) = 3;
    L(sub2indLP(sub)) = 1;
    sub(i) = 2;
end

end

function ind = sub2indLP(sub)
ind = 0;
for j = numel(sub):-1:1
    ind = ind * 3 + sub(j) - 1;
end
ind = ind + 1;
end