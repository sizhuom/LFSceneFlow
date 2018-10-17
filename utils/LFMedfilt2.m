function [ B ] = LFMedfilt2( A, msize )
%LFMEDFILT2 Apply median filtering on each 2-D slice of the light field

B = zeros(size(A));

for i = 1:size(A, 1)
    for j = 1:size(A, 2)
        for ch = 1:size(A, 5)
            B(i,j,:,:,ch) = medfilt2(squeeze(A(i,j,:,:,ch)), msize, 'symmetric');
        end
    end
end


end

