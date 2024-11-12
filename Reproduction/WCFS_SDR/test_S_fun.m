clc;
clear;
% % Test Code
N = 128;
M = 4;
threshold = 10;
fs = 128*1e6;
fc = 0e6;
B = 40e6;
T = 1e-6;
signal = generator_LFM(fs,fc,B,T); % 参考信号
% signal = ones(16,1);
% N = 16;
alpha = -(N-M)/2:(N-M)/2;
alpha = -N/2:N/2;
f = 1:N;
FN = dftmtx(N); 
Fshift = zeros(N);
Fshift(1:N/2,N/2+1:end) = eye(N/2);
Fshift(N/2+1:end,1:N/2) = eye(N/2);
FN = Fshift*FN;
S1 = zeros(size(alpha,2),size(f,2));
F = fs/(2*N);               % precompute constants -  F = fs/(2*N);     
G = fs/N; 
for k = 1:N
    k1 = 1:N;
    f(k,k1) = F*(k+k1-1) - fs/2;
    alfa(k,k1) = G*(k-k1 + N-1) - fs;
end
i_1 = 0; 
for i_a = alpha
    i_1 = i_1 + 1;
    i_2 = 0;
    for i_f = f
        i_2 = i_2 + 1;
        temp_taf = zeros(N,N);
        for i_t1 = 1:N
            for i_t2 = 1:N
                if i_t1 - i_t2 == i_a && i_t2 >=i_2 && i_t2 < i_2 + M
                    temp_taf(i_t1,i_t2) = 1;
                end
            end
        end

        S1(i_1, i_2) = (FN*signal)'*temp_taf*(FN*signal); 
    end
end

S3 = zeros(N,N);
% FNr
FN = dftmtx(N); 
Fshift = zeros(N);
Fshift(1:N/2,N/2+1:end) = eye(N/2);
Fshift(N/2+1:end,1:N/2) = eye(N/2);
FN = Fshift*FN;
Taf = zeros(N);
temp_m = -M/2+1:M/2; 
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
        S3(temp_1,temp_2) = (FN*signal)'*temp_taf*(FN*signal); 
    end
end
figure
plot(diag(abs(fliplr(S3)))); 


figure
mesh(abs(S1))
figure
mesh(f,alfa,abs(S3))
S2 = Analysis_CS_DFSM(fs, signal, 0, 4);
figure
plot(real(diag(S3)));
figure
mesh(abs(S3))


function A = setDiagonal(A, j_index, M, offset)
    % 获取斜对角线的长度
    len = length(diag(A, offset));
    
    % 创建一个值为value的向量
    values = zeros(len, 1);
    values(j_index:j_index+M-1) = 1;
    % 将向量设置为斜对角线的值
    A = spdiags(values, offset, A);
end