function [ str ] = struct2str( s )
%STRUCT2STR Simple struct to string function. No recursion.

fields = fieldnames(s);
str = '';
for i = 1:numel(fields)
    val = s.(fields{i});
    str = [str fields{i} '=' mat2str(val) ' '];
end

end

