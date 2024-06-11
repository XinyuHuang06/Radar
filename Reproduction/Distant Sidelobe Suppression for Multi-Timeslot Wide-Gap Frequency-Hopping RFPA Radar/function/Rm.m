% Example: Rm = Rm(MW, N1, N2, m, fnSeq, fpass)
% :param :
% :return :
% detailed description: 生成旁瓣能量计算矩阵
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 06/06/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function Rm = Rm(MW, N1, N2, m, fnSeq, fpass)
    Rm = zeros(MW);
    Nm = [m-N2:m-1,m+1:m+N1];
    fm = fnSeq(m);
    for i_1 = 1:MW
        for i_2 = 1:MW
            for i_Nm = Nm
                fn = fnSeq(i_Nm);
                Rm(i_1, i_2) = Rm(i_1, i_2) + ( 2*fpass*exp(1j*2*pi*(fn-fm)*(i_1-i_2)*TC)*sinc(2*fpass*(i_1-i_2)*TC) - ...
                                                fpass*exp(1j*2*pi*(fn-fm)*(i_1-i_2+1)*TC)*sinc(2*fpass*(i_1-i_2+1)*TC)-...
                                                fpass*exp(1j*2*pi*(fn-fm)*(i_1-i_2-1)*TC)*sinc(2*fpass*(i_1-i_2-1)*TC))...
                                                /(2*pi^2*(fn-fm)^2);
            end
        end
    end
end