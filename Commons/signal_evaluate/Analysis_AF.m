function [ambi,f_grid_plot, t_grid_plot] = Analysis_AF(varargin)
% Example: [ambi,f_grid_plot, t_grid_plot] = AF_Analysis(signal_1,signal_2,fs,'bool_mesh',1)
% :param signal_1: A N*1 column vector 
% :param signal_2: A N*1 column vector
% :param fs: the sampling rate
% :param bool_mesh: bool parameter, if it is true, the figure will be plot using the fun-mesh.
% :return t_grid_plot: the time-line sequence
% :return f_grid_plot: the frequenct-line sequence
% :return ambi: the normalized ambiguity function 
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 12/10/2023.
% Copyright (C) 2023 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
    in_par = inputParser;% Initialization
    addOptional(in_par, 'signal_1', 0);
    addOptional(in_par, 'signal_2', 0);
    addOptional(in_par, 'fs', 0);
    addParameter(in_par, 'bool_mesh' , 0);
    parse(in_par, varargin{:});  
    signal_1 = in_par.Results.signal_1;
    signal_2 = in_par.Results.signal_2;
    fs = in_par.Results.fs;
    bool_mesh = in_par.Results.bool_mesh;
    % Processing as column vectors
    if ~iscolumn(signal_1) 
        signal_1 = transpose(signal_1);
    end
    if ~iscolumn(signal_2) 
        signal_2 = transpose(signal_2);
    end
    % the number of symbols of signal
    N = numel(signal_1);
    % Time (frequency) normalized w.r.t. T (1/T)
    T = 1; 
    f_max = 0.5*fs;
    % the base sampling interval is T/N for time and 1/T for frequency
    osr = 1; % the over sampleing rate of the ambiguous function
    % Time grid
    t_grid_size = (T/N) / osr;
    Nt = N * osr;
    t_grid = (-Nt+1:Nt-1) * t_grid_size;
    % Fre grid
    % f_grid_size = (1/T) / osr;
    % Nf = ceil(f_max / f_grid_size);
    f_grid_size = f_max/1000;
    Nf = 1000;
    f_grid = (-Nf+1:Nf-1) * f_grid_size;
    % over sample signal
    m = reshape(ones(osr,1)*signal_1.',[1 Nt]);
    n = reshape(ones(osr,1)*signal_2.',[1 Nt]);
    % Calculate the AmbiFun
    ambi = zeros(2*Nt-1, 2*Nf-1);
    for iter_Df = 1:numel(f_grid)
        ambi(:,iter_Df) = xcorr(n,m.*exp(1j*2*pi*f_grid(iter_Df)*(0:Nt-1)*t_grid_size*N/fs));
    end
    ambi = abs(ambi);
    ambi = ambi / max(ambi(:));
    t_grid_plot = [-T, t_grid, T];
    f_grid_plot = [-f_max, f_grid, f_max];
    ambi = [zeros(1, 2*Nf-1); ambi; zeros(1, 2*Nf-1)];
    ambi = [zeros(2*Nt+1, 1), ambi, zeros(2*Nt+1, 1)];
    % Plot 3D AmbiFun
    if bool_mesh
        mesh(f_grid_plot, t_grid_plot, ambi);
        view(-30,30);
        xlabel('Doppler Shift/Hz');
        ylabel('Delay/\mus');
        zlabel('Normalized amplitude');
        SetDrawStyle;
    end
end