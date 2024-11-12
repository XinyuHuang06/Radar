% Example:
% :param :
% :return :
% detailed description:
% HEX2RGB: 实现颜色HEX值转化RGB
% 输入HEX的字符串'XXXXXX', 返回RGB值
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 29/02/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function RGB=HEX2RGB(HEX)
    H=['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','0']; % Hexadecimal tab
    RGB = zeros(3);
    for i=1:3
        y(1)=find(H==HEX(2*i-1)); % 在H中查找第2i-1个字符的位置     
        y(2)=find(H==HEX(2*i)); % 在H中查找第2i个字符的位置
        RGB(i)=y(1)*16+y(2); % 用公式计算RGB的值
    end
end

