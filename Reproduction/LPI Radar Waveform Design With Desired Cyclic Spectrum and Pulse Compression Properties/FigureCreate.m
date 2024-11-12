% Example: 
% :param : 
% :return : 
% detailed description: 
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 20/06/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
addpath("./function/");
RecordPath = ["output/01/Record.mat" "output/02/Record.mat" "output/03/Record.mat"];
% Fig 1
load(RecordPath(1));
DataRecordPack01 = DataRecordPack;
load(RecordPath(2));
DataRecordPack02 = DataRecordPack;
load(RecordPath(3));
DataRecordPack03 = DataRecordPack;
% 绘制不同参数下的对比图像
JsonPath = './data/initial_Parameter.json';
InitialParameter = loadjson(JsonPath);
fs = InitialParameter.signal.fs;
N = 256;
M = 4;
xr01 = DataRecordPack01.ParaRecord.xr{end};
xr02 = DataRecordPack02.ParaRecord.xr{end};
xr03 = DataRecordPack03.ParaRecord.xr{end};
cr = DataRecordPack01.ConstantPara.cr;
x01 = xr01(1:N) + 1j*xr01(N+1:end);
x02 = xr02(1:N) + 1j*xr02(N+1:end);
x03 = xr03(1:N) + 1j*xr03(N+1:end);
c = cr(1:N) + 1j*cr(N+1:end);
S_cr = Analysis_CS_DFSM(fs,cr,fs/N,M,'bool_draw',0);
S_cr = abs(CF_diag(c, N, M));
S_xr01 = abs(CF_diag(x01, N, M));
S_xr02 = abs(CF_diag(x02, N, M));
S_xr03 = abs(CF_diag(x03, N, M));
normal_value = max(max(S_cr));
profile_cr = diag(fliplr(S_cr/normal_value));
profile_xr01 = diag(fliplr(S_cr/normal_value));
profile_xr02 = diag(fliplr(S_cr/normal_value));
profile_xr03 = diag(fliplr(S_cr/normal_value));
fig1 = figure;
plot(profile_cr);hold on;
plot(profile_xr01);hold on;
plot(profile_xr02);hold on;
plot(profile_xr03);hold off;
title('\alpha = 0'); legend('cr','xr(\vartheta =0.1)','xr(\vartheta =0.2)','xr(\vartheta =0.3)');
SetDrawStyle;
exportgraphics(fig1, "./output/fig2.pdf");
% Fig 2