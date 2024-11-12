function varargout = PlotCompareFIgure(varargin)
    OutPath = varargin{1};
    load(strcat(OutPath(1),"/Record.mat"));
    cr = DataRecordPack.ConstantPara.cr;
    fs = DataRecordPack.ConstantPara.signal.fs;
    N  = DataRecordPack.ConstantPara.signal.N;
    M  = DataRecordPack.ConstantPara.signal.M;
    S_cr = CF_diag(cr(1:N)+1j*cr(N+1:end), N, M);
    normal_value = max(diag(abs(fliplr(S_cr))));
    profile_cr = diag(abs(fliplr(S_cr)));
    plot(profile_cr/normal_value);hold on;
    for i = 1:3
        load(strcat(OutPath(i),"/Record.mat"));
        if DataRecordPack.StopFlag
            xr = DataRecordPack.ParaRecord.xr{DataRecordPack.StopNum};
        else
            xr = DataRecordPack.ParaRecord.xr{end};
        end
        % 进行参数对比
        
        max(abs(xr-cr))
        out = CF_diag(xr(1:N)+1j*xr(N+1:end), N, M);
        out = abs(out);
        out = out/max(max(out));
        profile_cr = diag(abs(fliplr(out)));
        plot(profile_cr); 
    end
    title('\alpha = 0'); legend(["cr","xr(\delta=0.1)","xr(\delta=0.2)","xr(\delta=0.3)"]);
    SetDrawStyle;
end