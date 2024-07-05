function S = GenerateMatrixS(N, M, threshold)
    % FNr
    FN = dftmtx(N); 
    Fshift = zeros(N);
    Fshift(1:N/2,N/2+1:end) = eye(N/2);
    Fshift(N/2+1:end,1:N/2) = eye(N/2);
    FN = Fshift*FN;
    % 计算循环谱矩阵
    Taf = zeros(N);
    temp_m = -M/2+1:M/2; 
    for temp_1 = 1:N
        for temp_2 = 1:N
            temp_B = max(1-temp_1, 1-temp_2);
            temp_A = min (N-temp_1, N-temp_2);
            temp_n = temp_m((temp_m<=temp_A) & (temp_m>=temp_B));
            temp_taf = zeros(N,N);
            if isempty(temp_n)
                temp_taf = zeros(N,N);
            else
                temp_p = temp_1 + temp_n;
                temp_q = temp_2 + temp_n;
                temp_taf = zeros(N,N);
                temp_taf(temp_p(1):temp_p(end),temp_q(1):temp_q(end)) = eye(length(temp_q));
            end
            if abs(N + 1 - temp_1-temp_2) < threshold
                Taf = temp_taf*0.1 + Taf;
            else
                Taf = temp_taf*10 + Taf;
            end
        end
    end
    S = FN'*Taf*FN;
end