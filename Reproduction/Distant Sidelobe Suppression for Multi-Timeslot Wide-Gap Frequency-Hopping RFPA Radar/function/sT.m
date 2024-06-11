% Example: value = sT(t, tm, fm, phi0, TW, 'type', 'phasecode', 'TC', TC, 'z_m', z_m) 
%          value = sT(t, tm, fm, phi0, TW) 
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
function value = sT(varargin)
    in_par = inputParser;
    addOptional(in_par,'t', 0);
    addOptional(in_par, 'tm', 0);
    addOptional(in_par, 'fm', 0);
    addOptional(in_par, 'phi0', 0);
    addOptional(in_par, 'TW', 0);
    addParameter(in_par, 'type', 'normal');
    addParameter(in_par, 'TC', 0);
    addParameter(in_par, 'z_m', 0);
    parse(in_par, varargin{:}); 
    t = in_par.Results.t;
    tm = in_par.Results.tm;
    fm = in_par.Results.fm;
    phi0 = in_par.Results.phi0;
    TW = in_par.Results.TW;
    type = in_par.Results.type;
    TC = in_par.Results.TC;
    z_m = in_par.Results.z_m;
    if type == "normal"
        value = g(t-tm, TW).*exp(1j*(2*pi*fm*t + phi0));
    elseif type == "phasecode"
        value = g_pc(t, TW, TC, z_m).*exp(1j*(2*pi*fm*t + phi0));
    end
end

function value = g(t, TW)
    value = u(t/TW);
end

function value = g_pc(t, TW, TC, z_m)
    MW = TW/TC;
    value =  zeros(size(t));
    for i_MW = 1:MW
        value = z_m(i_MW)*u( (t-i_MW*TC)/TC );
    end
end

function value = u(t)
    value = zeros(size(t));
    value(t>=0 & t<=1) = 1;
end