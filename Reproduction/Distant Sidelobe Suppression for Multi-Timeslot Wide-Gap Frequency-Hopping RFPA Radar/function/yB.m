% Example: value = yB(M, max_l, max_k, xlk_Seq, n,p,TS,Delta_r,Delta_v,c,phi0,TW,tnSeq,fnSeq)
% :param :
% :return :
% detailed description: 返回添加了噪声之后的回波信号
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 06/06/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function value = yB(M, TargetMatrix, n,p,TS,Delta_r,Delta_v,c,phi0,TW,tnSeq,fnSeq)
    value = zeros(size(p));
    noise = rand(size(p)) + 1j*randn(size(p));
    Target_num = size(TargetMatrix, 2); % 第二维度为目标数目
    for i_M = 0:M-1
        for i_tar = 1:Target_num
            value = value + TargetMatrix(3,i_tar)*smlk(i_M,TargetMatrix(1,i_tar),TargetMatrix(2,i_tar),n,p,TS,Delta_r,Delta_v,c,phi0,TW,tnSeq,fnSeq);
        end
    end
    value = value + noise;
end