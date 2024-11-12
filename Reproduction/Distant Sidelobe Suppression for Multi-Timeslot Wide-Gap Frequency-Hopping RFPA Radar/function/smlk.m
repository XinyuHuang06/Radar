% Example: value = smlk(m, n, p, num_tar, Parameter);
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
function value = smlk(m, n, p, num_tar, Parameter)
    % Unpacket The parameter
    Delta_r = Parameter.Delta_r;
    Delta_v = Parameter.Delta_v;
    c = Parameter.c;
    phi0 = Parameter.phi0;
    tnSeq = Parameter.tnSeq;
    fnSeq = Parameter.fnSeq;
    TS = Parameter.TS;
    % Gneralize the signal
    l = Parameter.TargetMatrix(1,num_tar);
    k = Parameter.TargetMatrix(2,num_tar);
    p = p(:);
    tn = tnSeq(n);
    fn = fnSeq(n);
    tnp = tn + p*TS;
    if Parameter.flag_PC
        value = sT(tnp - 2*(l*Delta_r-k*Delta_v*tnp)/c, m, Parameter).*exp(-1j*(2*pi*fn*tnp+phi0+2));
    else
        value = sT(tnp - 2*(l*Delta_r-k*Delta_v*tnp)/c, m, Parameter).*exp(-1j*(2*pi*fn*tnp+phi0));
    end
    % value = sT(tnp - 2*(l*Delta_r-k*Delta_v*tnp)/c, m, Parameter);
%     value = sT( tnp - 2*(l*Delta_r-k*Delta_v*tnp)/c,  tm, fm, phi0, TW).*exp(-1j*(2*pi*fn*tnp+phi0));
end