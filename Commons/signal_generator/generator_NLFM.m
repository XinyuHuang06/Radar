function [signal_NLFM, t_grid] = generator_NLFM(varargin)
% Example: 
% :param fs:
% :param fc:
% :param 
% :return signal_NLFM:
% :return t:
% ������λ����ԭ�����NLFM�ź�
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 20/10/2023.
% Copyright (C) 2023 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
    in_par = inputParser;
    addOptional(in_par, 'fs', 0); % ����Ƶ��
    addOptional(in_par, 'fc', 0); % �ز�Ƶ��
    addOptional(in_par, 'B', 0); % ����
    addOptional(in_par, 'N', 0); % �źŵ���
    addOptional(in_par, 'fun', 0); % Ƶ���źź���
    parse(in_par, varargin{:});
    fs = in_par.Results.fs;
    fc = in_par.Results.fc;
    B = in_par.Results.B;
    N = in_par.Results.N;
    fun = in_par.Results.fun;

    T0 = N/fs;
    t_grid = (0:1/fs:(N-1)/fs)';
    fun = fun(:)/max(fun); % ת��Ϊ����������һ��
    fun = fun*B;
    m_fun = tril(ones(N))*fun; % ����
    m_fun = fun;
    carrier = exp(1j*2*pi*fc*t_grid);
    signal_NLFM = carrier.*exp(1j*2*pi*m_fun.*t_grid);
    
end