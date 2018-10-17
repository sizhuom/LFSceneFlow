function AbsRectH = IntrinRelToAbs( RelRectH )

AbsRectH = RelRectH;
AbsRectH(3:4,:) = RelRectH(3:4,:) + RelRectH(1:2,:);