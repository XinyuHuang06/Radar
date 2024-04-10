classdef forWaitbar < handle
    properties (Access = private)
        solid_square = '*'; % 已完成
        hollow_square = '-'; % 未完成
        N; % 总循环数
        temp_num; % 当前循环数
        square_nums = 40; % 打印方块数
        threshold; % 打印阈值
        Count_p_num; % 打印计数器
        s_num; % 方块计数器1
        h_num; % 方块计数器2
    end

    methods
        function obj = forWaitbar(N)
        % 初始化
            obj.N = N;
            obj.temp_num = 0;
            obj.Count_p_num = 0;
            obj.s_num = 0;
            obj.h_num = obj.square_nums;
            obj.threshold = fix(N/obj.square_nums);
        end
        function show_bar(obj)
            obj.Count_p_num = obj.Count_p_num + 1;
            obj.temp_num = obj.temp_num + 1;
            if obj.Count_p_num >= obj.threshold || obj.temp_num == obj.N || obj.temp_num == 1
                if obj.s_num ~= obj.square_nums - 1 || obj.temp_num == obj.N
                    obj.Count_p_num = 0;
                    obj.s_num = obj.s_num + 1;
                    obj.h_num = obj.h_num - 1;
                    obj.draw_bar;
                end
            end

        end
    end
    methods (Access = private)
        function draw_bar(obj)
            back = repmat('\b', [1,obj.square_nums+2+4]);
            fprintf(back);
            fprintf('<');
            if (obj.s_num > 0)
                for i = 1:obj.s_num
                    fprintf(obj.solid_square);
                end
            end
            if (obj.h_num > 0)
                for i = 1:obj.h_num
                    fprintf(obj.hollow_square);
                end
            end
            fprintf('>');
            precent = num2str(round((obj.temp_num/obj.N)*100));
            if length(precent) < 3
                precent = [repmat(' ',[1,3-length(precent)]),precent];
            end
            fprintf([precent,'%%']);
            if obj.temp_num == obj.N
                fprintf('\n');
            end
        end
    end
end