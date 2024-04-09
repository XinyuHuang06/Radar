
Num = 1000
N = 512;
M = 32;
fs = 10e6;
df = fs/10;
Sout = zeros(N,N);

for i = 1:Num
    signal = exp(1j*2*pi*rand(N,1))/sqrt(N);  
    [S, ~, ~] = Analysis_CS_DFSM(signal,fs, df, M);
    Sout = Sout + S/Num;

end

% 
% figure 
% Analysis_CS_FAM(fs,xr(1:N),fs/N,M,'bool_draw',1);
% exportgraphics(gcf, './output_files/xr_CS.pdf','ContentType', 'vector');
% 
% figure 
% Analysis_CS_DFSM(cr(1:N),fs,fs/N,M);
% exportgraphics(gcf, './output_files/cr_CS.pdf','ContentType', 'vector');