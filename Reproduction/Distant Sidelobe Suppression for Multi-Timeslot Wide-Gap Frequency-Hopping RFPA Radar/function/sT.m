% Example: value = sT(t, tm, Parameter)
% :param t: The current time, value or vector
% :param m: The No. of signal
% :param Parameter: The packet of parameter produced to the function of Initialization_Paramter.m
% :return value: The complex signal value, value or vector (the same of input parameter t)
% detailed description: To generalize the transmit siganl.
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 06/06/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function value = sT(t, m, Parameter)
    t = t(:);
    phi0 = Parameter.phi0;
    TW = Parameter.TW;
    tm = Parameter.tnSeq(m);
    fm = Parameter.fnSeq(m);
    if Parameter.flag_PC
        TC = Parameter.PhaseCodeParameter.TC;
        z_m = Parameter.PhaseCodeParameter.PhaseCode;
        value = g_pc(t-tm, TW, TC, z_m).*exp(1j*(2*pi*fm*t + phi0));  
    else
        value = g(t-tm, TW).*exp(1j*(2*pi*fm*t + phi0));
    end
end

function value = g(t, TW)
    value = u(t/TW);
end

function value = g_pc(t, TW, TC, z_m)
    MW = TW/TC;
    value = zeros(size(t));
    % 1
    index = floor(t(t>=0 & t<=TW)*MW/TW);
    if ~isempty(index)
        value(t>=0 & t<=TW) = z_m(index+1).*u((t(t>=0 & t<=TW)-index*TC)/TC);
    end
    % 2
    value(t<0 & t>TW) = 0;
end

function value = u(t)
    value = zeros(size(t));
    value(t>=0 & t<=1) = 1;
end