function PlotAndExport(DataSetPackets, DataRecordPack, InitialParameter, OutFolderPath)
    % Unpacket DataSet
    xr = DataSetPackets.xr;
    br = DataSetPackets.br;
    cr = DataSetPackets.cr;
    N = DataSetPackets.N;
    % Unpacket InitialParameter
    fs = InitialParameter.signal.fs;
    M = InitialParameter.signal.M;
    % Unpacket DataRecordPack
    TarS = DataRecordPack.TarRecord;
    if InitialParameter.Flag.Plot
        if exist(OutFolderPath,"dir") == 0
            mkdir(OutFolderPath);
        end
        % % The Sidelobe figure
        fig2 = figure;
        Analysis_Sidelobe(xr(1:N),xr(1:N),'bool_draw',1);hold on;
        Analysis_Sidelobe(cr(1:N),cr(1:N),'bool_draw',1);hold off;
        legend('xr','cr');
        SetDrawStyle;
        % siganl time-domain figure
        fig3 = figure;
        tiledlayout(2,2);
        nexttile;plot(xr);hold on;plot(cr);hold off;legend('xr','cr'); title('xr and cr');
        nexttile;plot(xr);hold on;plot(br);hold off;legend('xr','br'); title('xr and br');
        nexttile;plot(br); title('br');
        nexttile;plot(cr); title('cr');
        SetDrawStyle; 
        % CS DFSM
        fig4 = figure;
        tiledlayout(1,2);
        nexttile;[S_xr, ~, ~] = Analysis_CS_DFSM(fs,xr(1:N),fs/N,M,'bool_draw',1);
        nexttile;[S_cr, ~, ~] = Analysis_CS_DFSM(fs,cr(1:N),fs/N,M,'bool_draw',1);
        SetDrawStyle;
        fig5 = figure;
        plot(diag(S_cr));hold on;plot(diag(S_xr));hold off;
        title('\alpha = 0'); legend('cr','xr');
        SetDrawStyle;

        if InitialParameter.Flag.ExportFigure
        % SetDrawStyle;
            exportgraphics(fig2, strcat('./',OutFolderPath,'/Fig2_Sidelobe.png'),'ContentType', 'image');
            exportgraphics(fig3, strcat('./',OutFolderPath,'/Fig3_TimeDomain.png'),'ContentType', 'image');
            exportgraphics(fig4, strcat('./',OutFolderPath,'/Fig4_.png'),'ContentType', 'image');
            exportgraphics(fig5, strcat('./',OutFolderPath,'/Fig5_CS_x.png'),'ContentType', 'image',"Resolution",300);
        end
    end

end