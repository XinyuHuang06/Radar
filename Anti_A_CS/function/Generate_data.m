function [Data] = Generate_data(N, M, flag_Sparse, threshhold)
    % FNr
    FN = dftmtx(N); 
    Fshift = zeros(N);
    Fshift(1:N/2,N/2+1:end) = eye(N/2);
    Fshift(N/2+1:end,1:N/2) = eye(N/2);
    FN = Fshift*FN;
    FNr = complex2real(FN);
    % omega_alpha_1
    omega_alpha_1 = ones((N+1)*N/2,1);
    if flag_Sparse
        chi_matrix = cell(N,1);
        for i_temp = 1:N
            matrix_temp = zeros(N,N);
            matrix_temp(i_temp,i_temp) = 1;
            chi_matrix{i_temp} = sparse(complex2real(matrix_temp));
        end
        % 循环谱计算矩阵 20240410 Modified   
        temp_m = -M/2+1:M/2;
        Taf_1 = cell((N+1)*N/2,1);
        Taf_2 = cell((N+1)*N/2,1);
        Taf = cell((N+1)*N/2,1);
        i_temp = 0;
        for temp_1 = 1:N
            for temp_2 = temp_1:N
                i_temp = i_temp + 1;
                temp_B = max(1-temp_1, 1-temp_2);
                temp_A = min (N-temp_1, N-temp_2);
                temp_n = temp_m((temp_m<=temp_A) & (temp_m>=temp_B));
                if isempty(temp_n)
                    temp_taf = zeros(N,N);
                else
                    temp_p = temp_1 + temp_n;
                    temp_q = temp_2 + temp_n;
                    temp_taf = zeros(N,N);
                    temp_taf(temp_p(1):temp_p(end),temp_q(1):temp_q(end)) = eye(length(temp_q));
                end
                temp_taf_1 = [real(temp_taf),-imag(temp_taf);imag(temp_taf),real(temp_taf)];
                temp_taf_2 = [imag(temp_taf),-real(temp_taf);real(temp_taf),imag(temp_taf)];
                Taf_1{i_temp} = sparse(temp_taf_1);
                Taf_2{i_temp} = sparse(temp_taf_2);
                Taf{i_temp} = sparse(temp_taf);
                if abs(temp_1-temp_2) < threshhold
                    omega_alpha_1(i_temp) = 0.01;
                else
                    omega_alpha_1(i_temp) = 10;
                end
            end
        end    
    else
        chi_matrix = zeros(2*N,2*N,N); % 相似性约束矩阵
        for i_temp = 1:N
            matrix_temp = zeros(N,N);
            matrix_temp(i_temp,i_temp) = N^2;
            chi_matrix(:,:,i_temp) = complex2real(matrix_temp);
        end
        Taf_1 = zeros(2*N,2*N,(N+1)*N/2);
        Taf_2 = zeros(2*N,2*N,(N+1)*N/2);
        Taf = zeros(N,N,(N+1)*N/2);
        temp_m = -M/2+1:M/2;
        i_temp = 0;
        for temp_1 = 1:N
            for temp_2 = temp_1:N
                i_temp = i_temp + 1;
                temp_B = max(1-temp_1, 1-temp_2);
                temp_A = min (N-temp_1, N-temp_2);
                temp_n = temp_m((temp_m<=temp_A) & (temp_m>=temp_B));
                if isempty(temp_n)
                    temp_taf = zeros(N,N);
                else
                    temp_p = temp_1 + temp_n;
                    temp_q = temp_2 + temp_n;
                    temp_taf = zeros(N,N);
                    temp_taf(temp_p(1):temp_p(end),temp_q(1):temp_q(end)) = eye(length(temp_q));
                end
                temp_taf_1 = [real(temp_taf),-imag(temp_taf);imag(temp_taf),real(temp_taf)];
                temp_taf_2 = [imag(temp_taf),-real(temp_taf);real(temp_taf),imag(temp_taf)];
                Taf_1(:,:,i_temp) = sparse(temp_taf_1);
                Taf_2(:,:,i_temp) = sparse(temp_taf_2);
                Taf(:,:,i_temp) = temp_taf;
                if abs(temp_1-temp_2) < threshhold
                    omega_alpha_1(i_temp) = 0.01;
                else
                    omega_alpha_1(i_temp) = 10;
                end
            end
        end  
    end
    % Intialize Data Struct
    Data = struct('chi_matrix',{chi_matrix},'Taf_1',{Taf_1},'Taf_2',{Taf_2},'FNr',{FNr},'omega_alpha',{omega_alpha_1},'Taf',{Taf},'sparse',{flag_Sparse});
    Data.N = N;
    Data.FN = FN;
end