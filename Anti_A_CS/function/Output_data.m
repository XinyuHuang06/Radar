function [] = Output_data(filename, TarS)
    TarSfields = fieldnames(TarS);
    % 如果文件已经存在，则删除
    if exist(filename,"file")
        delete(filename);
    end
    letters = 'ABCDEFGHYJKLMNOPQRSTUVWXYZ';
    for i_Tar = 1:length(TarSfields)
        key_t = TarSfields{i_Tar};
        value = TarS.(key_t);
        writematrix(key_t,filename,'Sheet',1,'Range',[letters(i_Tar),'1']);
        writematrix(zeros(),filename,'Sheet',1,'Range',[letters(i_Tar),'2']);
        writematrix(value,filename,'Sheet',1,'Range',[letters(i_Tar),'2']);
    end
end