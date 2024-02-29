function rho = PAPR(x)
% Example:PAPR(signal)
% :param x: the real vector signal
% :return : the peak-average-power-ratio of signal x
% detailed description:
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 09/11/2023.
% Copyright (C) 2023 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
    rho = (max(abs(x)))^2 / ((norm(x))^2 / length(x));
end
