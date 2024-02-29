function [out,PSL,ISL] = Analysis_Sidelobe(varargin)
    % Example:
    % :param :
    % :return :
    % detailed description:
    %------------------------------------------------------------------------------
    % Created by: Xinyu Huang.
    % On: 14/11/2023.
    % Copyright (C) 2023 Xinyu Huang (learning_huang@163.com).
    % All Rights Reserved.
    % Unauthorized copying of this file, via any medium is strictly prohibited.
    % Proprietary and confidential.
    %------------------------------------------------------------------------------
        in_par = inputParser;
        addOptional(in_par, 'signal_1', 0);
        addOptional(in_par, 'signal_2', 0);
        addParameter(in_par, 'bool_draw', 1);
        parse(in_par,varargin{:});
        x1 = in_par.Results.signal_1;
        x2 = in_par.Results.signal_2;
        bool_draw = in_par.Results.bool_draw;
        x1 = x1(:);
        x2 = x2(:);
        N1 = length(x1);
        N2 = length(x2);
        N = max(N1, N2);
        if N1>N2
            x2 = [x2;zeros(N1-N2,1)];
        elseif N2>N1
            x1 = [x1;zeros(N2-N1,1)];
        end
        out = xcorr(x1,x2,"normalized");
        out = abs(out);
        out = out/max(out);

        PSL = max((out(N+1:end)));
        ISL = sum((out(N+1:end)).^2);

        out = 10*log10(out); % dB
        if bool_draw
            n_grid_plot = -N+1:1:N-1;
            plot(n_grid_plot, out);
        end
    end