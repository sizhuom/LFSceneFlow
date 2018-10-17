function [ LFc ] = colorCorrect( LF, lfrFile )
%COLORCORRECT Color-correct light fields using LFColourCorrect

lfp = LFReadLFP(lfrFile);
LFMetadata = lfp.Metadata;
ColMatrix = reshape(LFMetadata.image.color.ccm,3,3);
ColBalance = [...
    LFMetadata.image.color.whiteBalanceGain.r, ...
    LFMetadata.image.color.whiteBalanceGain.gb, ...
    LFMetadata.image.color.whiteBalanceGain.b ];
Gamma = 1;
LFc = LFColourCorrect(LF,ColMatrix,ColBalance,Gamma);

end

