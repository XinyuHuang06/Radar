addpath('../Commons/signal_generator');

T = 100*1e-6;
B = 3*1e6;
fs = 4*B;
fc = 0e6;

[signal_LFM, t] = generator_LFM(fs,fc,B,T,'type',1);

fig1 = figure;
% plot(real(signal_LFM));
spectrogram(signal_LFM, 128, 127, 128, fs, 'yaxis');

exportgraphics(fig1, './output/LFM_simulation.jpg', 'Resolution', 300);

% [out_freseq,out_FreSpec_Ampli,out_FreSpec_Phase] = Analysis_FS(signal_LFM, fs, 1);
exportgraphics(gcf, './output/LFM_simulation2.jpg', 'Resolution', 300);