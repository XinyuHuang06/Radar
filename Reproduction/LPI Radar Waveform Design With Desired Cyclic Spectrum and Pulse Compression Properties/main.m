profile on
addpath("./function/");
addAllSubfolders("../../Commons/");
DataPath = ["output/parameter_01.mat" "output/parameter_02.mat" "output/parameter_03.mat"];
OutPath = ["output/01" "output/02" "output/03"];
%% Initial
Initial
%% Running
for i = 1:3
    WCFS_ADMM(DataPath(i), OutPath(i));
end
fig1 = figure;
PlotCompareFIgure(OutPath);
% exportgraphics(fig1, "./output/fig2.pdf");
profile viewer


