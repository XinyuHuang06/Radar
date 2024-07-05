addpath("./function/");

N = 3; % 对比源数量
DataPath = ["data/parameter_01.mat" "data/parameter_02.mat" "data/parameter_03.mat"];
OutPath = ["output/01" "output/02" "output/03"];

Initial
for i = 1:3
    WCFS_ADMM(DataPath(i), OutPath(i));
end


% % 绘制结果图像

load(strcat(OutPath(1),"/Record.mat"));
cr = DataRecordPack.ConstantPara.cr;
fs = DataRecordPack.ConstantPara.signal.fs;
N  = DataRecordPack.ConstantPara.signal.N;
M  = DataRecordPack.ConstantPara.signal.M;
S_cr = CF_diag(cr(1:N)+1j*cr(N+1:end), N, M);
fig1 = figure;
plot(diag(abs(fliplr(S_cr))));hold on;
for i = 1:3
    load(strcat(OutPath(i),"/Record.mat"));
    if DataRecordPack.StopFlag
        xr = DataRecordPack.ParaRecord.xr{DataRecordPack.StopNum};
    else
        xr = DataRecordPack.ParaRecord.xr{end};
    end
    max((xr-cr))
    out = CF_diag(xr(1:N)+1j*xr(N+1:end), N, M);
    plot(diag(abs(fliplr(out)))); 
end
title('\alpha = 0'); legend(["cr","xr(\delta=0.1)","xr(\delta=0.2)","xr(\delta=0.3)"]);
exportgraphics(fig1, "./output/fig24.pdf");
SetDrawStyle;
rmpath("function/");