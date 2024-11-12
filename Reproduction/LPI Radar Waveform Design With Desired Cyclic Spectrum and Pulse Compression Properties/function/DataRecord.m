classdef DataRecord < handle
    properties (Access = public)
        TarRecord;
        ParaRecord;
        ConstantPara;
        % Other 
        MaxItersNum;
        StopFlag;
        StopNum;
    end
    methods
        function obj = DataRecord(Max_ItersNum)
            % The target function value record.
            obj.TarRecord.SUM = zeros(Max_ItersNum+1,1);
            obj.TarRecord.TargetFun = zeros(Max_ItersNum+1,1);
            obj.TarRecord.TargetLagrange_1 = zeros(Max_ItersNum+1,1);
            obj.TarRecord.TargetLagrange_2 = zeros(Max_ItersNum+1,1);
            obj.TarRecord.TargetPenalty_1 = zeros(Max_ItersNum+1,1);
            obj.TarRecord.TargetPenalty_2 = zeros(Max_ItersNum+1,1);

            % The parameter data record.
            obj.ParaRecord.xr = cell(Max_ItersNum+1,1);
            obj.ParaRecord.rho_0 = cell(Max_ItersNum+1,1);
            obj.ParaRecord.rho_1 = cell(Max_ItersNum+1,1);
            obj.ParaRecord.lambda_0 = cell(Max_ItersNum+1,1);
            obj.ParaRecord.lambda_1 = cell(Max_ItersNum+1,1);
            obj.ParaRecord.h = cell(Max_ItersNum+1,1);
            % Constant Parameter
            obj.ConstantPara.cr = 0;
            % Other Parameter
            obj.MaxItersNum = Max_ItersNum;
            obj.StopFlag = false;
        end
    end
    methods (Access = public)

        function UpdateTarRecord(obj, num, TarRecord)
            obj.TarRecord.SUM(num) = TarRecord.SUM;
            obj.TarRecord.TargetFun(num) = TarRecord.TargetFun;
            obj.TarRecord.TargetLagrange_1(num) = TarRecord.TargetLagrange_1;
            obj.TarRecord.TargetLagrange_2(num) = TarRecord.TargetLagrange_2;
            obj.TarRecord.TargetPenalty_1(num)= TarRecord.TargetPenalty_1;
            obj.TarRecord.TargetPenalty_2(num) = TarRecord.TargetPenalty_2;
        end

        function UpdateParaRecord(obj, num, value, str, varargin)
            if nargin > 4
                updatekey(obj, num, value, str)
                M = (nargin-4)/2;
                for i = 1:M
                    value = varargin{(2*(i-1)+1)};
                    str = varargin{2*i};
                    updatekey(obj, num, value, str)
                end
            else
                updatekey(obj, num, value, str)
            end
        end
        
        function InitialConstantParameter(obj, JsonPath)
            InitialParameter = loadjson(JsonPath);
            obj.ConstantPara.Max_ItersNum = InitialParameter.Max_ItersNum;
            obj.ConstantPara.Threshold = InitialParameter.Threshold;
            obj.ConstantPara.Flag = InitialParameter.Flag;
            obj.ConstantPara.signal = InitialParameter.signal;
            obj.ConstantPara.signal.N = round(InitialParameter.signal.T*InitialParameter.signal.fs);
            signalpara = obj.ConstantPara.signal;
            c = generator_LFM(signalpara.fs,signalpara.fc,signalpara.B,signalpara.T);
            cr = [real(c);imag(c)];
            obj.ConstantPara.cr = cr;
        end
    
        function flag = JudgeStop(obj, iter_num, threshold)
            if abs( obj.TarRecord.SUM(iter_num) - obj.TarRecord.SUM(iter_num-1) ) < threshold
                obj.StopFlag = true;
                obj.StopNum = iter_num;
                flag = true;
            else
                flag = false;
            end
        end
    end
    methods (Access = private)
        function updatekey(obj, num, value, str)
            if any( contains(fieldnames(obj.ParaRecord),str) )
                obj.ParaRecord.(str){num} = value;
            else
                fprintf('Unknown keywords!\n');
            end
        end
    end
end