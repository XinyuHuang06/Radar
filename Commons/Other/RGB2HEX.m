% Example: 
% :param :
% :return :
% detailed description:
% RGB2HEX : 实现颜色RGB值转化HEX
% 输入RGB三个数的数组[a,b,c], 返回HEX值
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 29/02/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function HEX=RGB2HEX(RGB)
    H=['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','0'];
    for i=1:3
        y(1)=floor(RGB(i)/16);
        y(2)=mod(RGB(i),16);
        HEX(2*i-1)=H(y(1));
        HEX(2*i)=H(y(2));
    end
end

