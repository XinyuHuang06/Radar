% Example:
% :param :
% :return :
% detailed description:
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 06/06/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function tnSeq = generate_tnSeq(M,TJ,TR)
    base_tn = ((1:M)*TR)'; % 平均PRI间隔
    tnSeq = base_tn + TJ*rand(M,1); % 偏移
end