
function [] = Plot_cruve(xr,br,cr,fs,N,M,TarS,flag_PlotAndExport)
    % i_m = length(TarS.Tar)-1;
    if flag_PlotAndExport 
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
        fig4 = figure;
        tiledlayout(1,2);
        nexttile;[S_xr, ~, ~] = Analysis_CS_DFSM(fs,xr(1:N),fs/N,M,'bool_draw',1);
        nexttile;[S_cr, ~, ~] = Analysis_CS_DFSM(fs,cr(1:N),fs/N,M,'bool_draw',1);
        SetDrawStyle;
        fig5 = figure;
        % plot(diag(S_cr(1:end,end:-1:1)));hold on;
        % plot(diag(S_xr(1:end,end:-1:1)));hold off;
        plot(diag(S_cr));hold on;
        plot(diag(S_xr));hold off;
        title('\alpha = 0')
        legend('cr','xr');
        SetDrawStyle;
        % fig7 = figure;
        % [ambi,~,~] = Analysis_AF(xr(1:N),xr(1:N),fs,'bool_mesh',1);SetDrawStyle;
        % fig8 = figure;
        % [~,~,~] = Analysis_FS(xr(1:N),fs,1);SetDrawStyle;
        % fig9 = figure;
        % [~,~,~] = Analysis_FS(cr(1:N),fs,1);SetDrawStyle;

    
        fig10 = figure;
        TarSfields = fieldnames(TarS);
        lenged_str = cell(length(TarSfields),1);
        for i_Tar = 1:length(TarSfields)
            key_t = TarSfields{i_Tar};
            value = TarS.(key_t);
            plot(0:length(value)-1,value); hold on;
            lenged_str{i_Tar} = key_t;
        end
        hold off;
        legend(lenged_str);
        SetDrawStyle;
        exportgraphics(fig2, './output_files/Fig2_Sidelobe.png','ContentType', 'image');
        exportgraphics(fig3, './output_files/Fig3_TimeDomain.png','ContentType', 'image');
        exportgraphics(fig4, './output_files/Fig4_.png','ContentType', 'image');
        exportgraphics(fig5, './output_files/Fig5_CS_x.png','ContentType', 'image',"Resolution",300);
        % exportgraphics(fig7, './output_files/Fig7_AF_x.png','ContentType', 'image');
        % exportgraphics(fig8, './output_files/Fig8_FS.png','ContentType', 'image');
        % exportgraphics(fig9, './output_files/Fig8_FS_2.png','ContentType', 'image');
        exportgraphics(fig10, './output_files/Fig10_Tar.png','ContentType', 'image');
    end

end