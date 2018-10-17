function [ LFgray ] = LFRGB2Gray( LFrgb )
%LFRGB2GRAY Convert a 4D RGB light field to grayscale
if (size(LFrgb,5) > 1)
    LFgray = 0.299 * LFrgb(:, :, :, :, 1) + 0.587 * LFrgb(:, :, :, :, 2) + ...
        0.114 * LFrgb(:, :, :, :, 3);
else
    LFgray = LFrgb;
end

end

