function [ fname ] = findNextTmp( workDir, pattern )
%FINDNEXTTMP Find next available filename for a temp file
% Assume '*', '%04d'

listing = dir(fullfile(workDir,pattern));
ind = find(pattern == '*');
ind = ind(1);
scanPat = [pattern(1:ind-1) '%04d' pattern(ind+1:end)];
max = 0;
for i = 1:length(listing)
    num = sscanf(listing(i).name, scanPat);
    if num > max
        max = num;
    end
end

fname = fullfile(workDir,sprintf(scanPat, max+1));

end

