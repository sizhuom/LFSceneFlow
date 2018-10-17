function [ str ] = cell2str( c )
%cell2STR Simple cell to string function. No recursion.

str = '';
for i = 1:numel(c)
    str = [str mat2str(c{i}) ' '];
end

end

