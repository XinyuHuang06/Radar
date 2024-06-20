addpath("function/");
DataPath = ["data/parameter_01.mat" "data/parameter_02.mat" "data/parameter_03.mat"];
OutPath = ["output/01" "output/02" "output/03"];
parfor i = 1:3
    WCFS_ADMM(DataPath(i), OutPath(i));
end
rmpath("function/");