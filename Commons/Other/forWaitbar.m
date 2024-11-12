classdef forWaitbar < handle
    properties (Access = private)
        solid_square = '*'; % The symbol of completion
        hollow_square = '-'; % The symbol of uncompletion
        N; % the numers of 
        temp_num; % 当前循环数
        square_nums = 40; % 打印方块数
        s_num; % solid_square acculation
        h_num; % hollow_square acculation
    end

    methods
        function obj = forWaitbar(N)
        % 初始化
            obj.N = N;
            obj.temp_num = 0;
            obj.s_num = 0;
            obj.h_num = obj.square_nums;
        end
        function show_bar(obj)
            obj.temp_num = obj.temp_num + 1;
            if obj.temp_num/obj.N >= obj.s_num/obj.square_nums
                if ( obj.s_num ~= obj.square_nums - 1 || obj.temp_num == obj.N ) && obj.temp_num <= obj.N
                    obj.s_num = rem(round(obj.temp_num*obj.square_nums/obj.N),obj.square_nums+1);
                    obj.h_num = obj.square_nums-obj.s_num;
                    obj.draw_bar;
                end
            end
        end
    end
    methods (Access = private)
        function draw_bar(obj)
            back = repmat('\b', [1,obj.square_nums+2+4]);
            fprintf(back);
            if (obj.s_num > 0)
                solid = repmat(obj.solid_square, [1,obj.s_num]);
            else
                solid = '';
            end
            if (obj.h_num > 0)
                hollow = repmat(obj.hollow_square, [1,obj.h_num]);
            else
                hollow = '';
            end
            precent = num2str(round((obj.temp_num/obj.N)*100));
            if length(precent) < 3
                precent = [repmat(' ',[1,3-length(precent)]),precent];
            end
            fprintf(['<',solid,hollow,'>',precent,'%%']);
            if obj.temp_num == obj.N
                fprintf('\n');
            end
        end
    end
end