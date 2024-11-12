function [] = print_log(TarS,s_xr_init,s_xr_post,s_cr,flag_PrintLog)
    if flag_PrintLog 
        try
            fid = fopen('./output_files/WCS_sparse.log','a+');
            date_str = char(datetime("now"));
            fprintf(fid, ['[',date_str,']','\n']);
            % % 输出波形参数信息
            print_signal_info(fid, s_xr_init);
            print_signal_info(fid, s_xr_post);
            print_signal_info(fid, s_cr);
            % % 输出相似度
            % xr_post and cr
            similiar_rate = (norm(s_cr.signal) - norm(s_xr_post.signal - s_cr.signal))/norm(s_cr.signal);
            fprintf(fid,['The similiar metric on xr and cr: %.2d\n'],similiar_rate);
            % % 输出迭代相关信息
            TarSfields = fieldnames(TarS);
            max_len = max(cellfun(@length, TarSfields));
            fprintf(fid,[repmat(' ', [1,max_len+1]),'Initial',repmat(' ', [1,6]),'end\n']);
            for i_Tar = 1:length(TarSfields)
                key_t = TarSfields{i_Tar};
                value = TarS.(key_t);
                if length(key_t) < max_len
                    add_str = repmat(' ', [1,max_len-length(key_t)]);
                else
                    add_str = '';
                end
                fprintf(fid,[key_t,add_str,': ','%.2d',repmat(' ', [1,5]),'%.2d.\n'],value(1),value(end));
            end
            fclose(fid);
        catch
            fprintf(fid,['print_log function error!']);
            fclose(fid);
        end
    end
end


function print_signal_info(fid, signal)
    fprintf(fid,[signal.s_name,'\n']);
    fprintf(fid,['PAPR: %.2f. ', ...
                'PSL: %.2f dB. ',...
                'ISL: %.2f dB. ',...
                'PSLR: %.2f dB. ',...
                'ISLR: %.2f dB.\n'...
                ],signal.PAPR,signal.PSL,signal.ISL,signal.PSLR,signal.ISLR);
end