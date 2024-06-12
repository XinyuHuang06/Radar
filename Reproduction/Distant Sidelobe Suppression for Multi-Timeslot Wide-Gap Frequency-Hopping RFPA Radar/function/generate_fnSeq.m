% Example: fnSeq = generate_fnSeq(Parameter);
% :param : Parameter结构体由Initialization_Parameter函数生成。
% :return :
% detailed description: 依照准则生成载频序列
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 06/06/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function fnSeq = generate_fnSeq(fC, G, B, N, M)
    fnSeq = zeros(M,1);
    for n=0:M-1
        while(1)
            f = fC + B*(rand-0.5);
            if n==0
                fnSeq(n+1) = f;
                break;
            else
                temp_1 = min(n,N);
                if min( abs(f-fnSeq(n-temp_1+1:n)) ) >= G
                    fnSeq(n+1) = f;
                    break;
                end
            end
        end
    end
end