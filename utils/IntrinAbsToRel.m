function RelRectH = IntrinAbsToRel( AbsRectH )

RelRectH = AbsRectH;
RelRectH(3:4,:) = AbsRectH(3:4,:) - AbsRectH(1:2,:);