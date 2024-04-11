num_Mator = 1000;
N = 512;
M = 4;
x =  exp(1j*2*pi*rand(N,1))/sqrt(N);
fs = 1e6;
iter_waitbar = forWaitbar(num_Mator);
Out = 0;
df = fs/100;
for i = 1:num_Mator
    [CS,f,alfa] = Analysis_CS_DFSM(fs,x,df,M,'bool_draw',0);
    Out = Out + CS/num_Mator;
    iter_waitbar.show_bar;
end
contour(alfa, f, CS); grid;
xlabel('Cycle frequency (Hz)'); ylabel('Frequency (Hz)');
title (['Frequency Smoothing SCD ', ', df = ', int2str(df),', N = ', int2str(N)]);
colormap(othercolor('PuBu7'));
colorbar;
SetDrawStyle;
exportgraphics(gcf, './output_files/ideal_CS.pdf','ContentType', 'vector');