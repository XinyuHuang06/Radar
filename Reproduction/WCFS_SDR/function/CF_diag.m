function out = CF_diag(signal, N, M)
    signal = (signal);
    FN = dftmtx(N); 
    Fshift = zeros(N);
    Fshift(1:N/2,N/2+1:end) = eye(N/2);
    Fshift(N/2+1:end,1:N/2) = eye(N/2);
    FN = Fshift*FN/sqrt(N);
    temp_m = -M/2+1:M/2; 
    out = zeros(N,N);
    for temp_1 = 1:N
        for temp_2 = 1:N
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
            out(temp_1,temp_2) = (FN*signal)'*temp_taf*(FN*signal); 
        end
    end
end