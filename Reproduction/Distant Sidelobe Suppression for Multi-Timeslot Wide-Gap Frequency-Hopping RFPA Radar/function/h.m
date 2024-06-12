% Example: value = h(Parameter);
% :param :
% :return :
% detailed description: 使用Kaiser窗生成一类FIR低通滤波器
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 06/06/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function value = h(fpass,fstop,fS,Rp,As)
    FS2 = fS/2; % 奈奎斯特频率 64 MHz
    f = [fpass fstop]/FS2;
    A = [1 0];
    % delta1 = (10^(Rp/20)-1)/(10^(Rp/20)+1);
    % delta2 = (1+delta1)*(10^(-As/20));
    % dev = [delta1 delta2];
    dev = [10^(-Rp/20) 10^(-As/20)];
    [N,Wn,beta,~] = kaiserord(f, A, dev);
    N = N + rem(N,2); % 确保滤波器阶数为奇数
    value = fir1(N, Wn, kaiser(N+1, beta));
end