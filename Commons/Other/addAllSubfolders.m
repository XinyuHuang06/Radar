function addAllSubfolders(folder)
    % 添加指定文件夹及其所有子文件夹到搜索路径中
    addpath(folder);
    subfolders = genpath(folder);
    addpath(subfolders);
end