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

% 利用已有的数据文件绘制性能对比图

% 数据文件路径

RecordPath = ["output/01/Record.mat" "output/02/Record.mat" "output/03/Record.mat"];


% 绘制对比图像

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



xr01 = DataRecordPack01.ParaRecord.xr{end};
xr02 = DataRecordPack02.ParaRecord.xr{end};
xr03 = DataRecordPack03.ParaRecord.xr{end};
x01 = xr01(1:N) + 1j*xr01(N+1:end);
x02 = xr02(1:N) + 1j*xr02(N+1:end);
x03 = xr03(1:N) + 1j*xr03(N+1:end);
figure(1);
S_xr01 = Analysis_CS_DFSM(fs,xr01,fs/N,M,'bool_draw',0);
S_xr02 = Analysis_CS_DFSM(fs,xr02,fs/N,M,'bool_draw',0);
S_xr03 = Analysis_CS_DFSM(fs,xr03,fs/N,M,'bool_draw',0);
S_cr = Analysis_CS_DFSM(fs,c,fs/N,M,'bool_draw',0);
plot(diag(S_cr));hold on;
plot(diag(S_xr01));hold on;
plot(diag(S_xr02));hold on;
plot(diag(S_xr03));hold off;
title('\alpha = 0'); legend('cr','xr(\vartheta =0.1)','xr(\vartheta =0.2)','xr(\vartheta =0.3)');
SetDrawStyle;
% Fig 2