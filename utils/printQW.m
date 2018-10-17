function printQW( arg1, arg2 )
%PRINTQW print the qw representation in millimeters and degrees

if nargin == 1
    fid = 1;
    qw = arg1;
else
    fid = arg1;
    qw = arg2;
end

for i = 1:size(qw, 2)
    fprintf(fid, '%10.06f ', qw(1:3, i)*1000);
    fprintf(fid, '%10.06f ', qw(4:6, i)/pi*180);
    fprintf(fid, '\n');
end


end

