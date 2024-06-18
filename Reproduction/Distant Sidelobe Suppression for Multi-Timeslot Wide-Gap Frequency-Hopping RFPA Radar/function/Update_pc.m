% Example: Update_pc(m, Rm, Parameter)
% :param : 
% :return : 
% detailed description: 
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 12/06/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function Update_pc(m, Rm,threshold, Parameter)
    % 更新第m个信号的相位编码
    Iter_MaxNum = 1000;
    MW = Parameter.PhaseCodeParameter.MW;
    pc_1 = Parameter.PhaseCodeParameter.PhaseCode(:,m);
    temp_pc = pc_1;
    % To update the phasecode by the Coordinate-Descent.
    dis  = 1e6;
    count = 0;
    while( dis > threshold )
        last_pc = temp_pc;
        for i = 1:MW
            temp_pc(i) = 0;
            e_i = zeros(MW,1);
            e_i(i) = 1;
            temp_code = -exp(1j*angle(e_i'*Rm*temp_pc));
            temp_pc(i) = temp_code;
        end
        count = count + 1;
        dis = norm(last_pc-temp_pc, 2);
        if count > Iter_MaxNum
            break;
        end
    end
    Parameter.PhaseCodeParameter.PhaseCode(:,m) = temp_code;
end