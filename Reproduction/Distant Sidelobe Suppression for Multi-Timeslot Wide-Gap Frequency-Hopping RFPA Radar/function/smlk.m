% Example: value = smlk(m,l,k,n,p,TS,Delta_r,Delta_v,c,phi0,TW,tnSeq,fnSeq)
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
function value = smlk(m,l,k,n,p,TS,Delta_r,Delta_v,c,phi0,TW,tnSeq,fnSeq)
    tn = tnSeq(n+1);
    tm = tnSeq(m+1);
    fn = fnSeq(n+1);
    fm = fnSeq(m+1);
    tnp = tn + p*TS;
    value = sT( tnp - 2*(l*Delta_r-k*Delta_v*tnp)/c,  tm, fm, phi0, TW).*exp(-1j*(2*pi*fn*tnp+phi0));
    % value = sT( tnp - 2*(l*Delta_r-k*Delta_v*tnp)/c,  tm, fm, phi0, TW).*exp(-1j*(2*pi*fn*tn+phi0));
end