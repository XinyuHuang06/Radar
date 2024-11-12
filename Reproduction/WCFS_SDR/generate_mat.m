addpath("function/")
fs = 128*1e6;
B = 20*1e6;
fc = 0;
T  = 2*1e-6;
N  = ceil(T*fs);
M  = 4;
t  = (0:1/fs:(N-1)/fs)';
f_ordinate = zeros(N,N);
alfa_ordinate = zeros(N,N);
for i_1 = 1:N
    i_2 = 1:N;
    f_ordinate(i_1,i_2) = (i_1 + i_2 -1)/(2*N) - 1/2;
    alfa_ordinate(i_1,i_2) = (i_1 - i_2 + N - 1)/N - 1;
end