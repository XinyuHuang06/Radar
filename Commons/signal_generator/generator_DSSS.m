function [signal_DSSS] = generator_DSSS(varargin)
% Example: DSSS_generator(signal,PNcode)
% :param signal: a N*1 column vector
% :param PNcode: a NSS*1 column vector, PNcode(i)-{-1,1}
% :return signal_DSSS: The signal after direct sequence spread spctrum
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 25/10/2023.
% Copyright (C) 2023 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
    in_par = inputParser; % Intialization
    addOptional(in_par, 'signal', 0);
    addOptional(in_par, 'PNcode', -1);
    parse(in_par, varargin{:});
    signal = in_par.Results.signal;
    PNcode = in_par.Results.PNcode;
    PNcode = PNcode(:);
    signal = signal(:);
    NSS = length(PNcode);
    N = length(signal);
    signal_DSSS = reshape(PNcode*transponse(signal),N*NSS,1);
end