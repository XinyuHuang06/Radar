classdef DataRecord < handle
    properties (Access = public)
        TarRecord;
        ParaRecord;
    end
    methods
        function obj = DataRecord(Max_ItersNum)
            obj.TarRecord.SUM = zeros(Max_ItersNum+1,1);
            obj.TarRecord.TargetFun = zeros(Max_ItersNum+1,1);
            obj.TarRecord.TargetLagrange_1 = zeros(Max_ItersNum+1,1);
            obj.TarRecord.TargetLagrange_2 = zeros(Max_ItersNum+1,1);
            obj.TarRecord.TargetPenalty_1 = zeros(Max_ItersNum+1,1);
            obj.TarRecord.TargetPenalty_2 = zeros(Max_ItersNum+1,1);


            obj.ParaRecord.xr = cell(Max_ItersNum+1,1);
            obj.ParaRecord.rho_0 = cell(Max_ItersNum+1,1);
            obj.ParaRecord.rho_1 = cell(Max_ItersNum+1,1);
            obj.ParaRecord.lambda_0 = cell(Max_ItersNum+1,1);
            obj.ParaRecord.lambda_1 = cell(Max_ItersNum+1,1);
        end
    end
    methods (Access = public)
        function UpdateTarRecord(obj, TarRecord, num)
            obj.TarRecord.SUM(num) = TarRecord.SUM;
            obj.TarRecord.TargetFun(num) = TarRecord.TargetFun;
            obj.TarRecord.TargetLagrange_1(num) = TarRecord.TargetLagrange_1;
            obj.TarRecord.TargetLagrange_2(num) = TarRecord.TargetLagrange_2;
            obj.TarRecord.TargetPenalty_1(num)= TarRecord.TargetPenalty_1;
            obj.TarRecord.TargetPenalty_2(num) = TarRecord.TargetPenalty_2;
        end
        function UpdateParaRecord(obj, value, str, varargin)
            if nargin > 3
                update_packet(obj, value, str);
                M = (nargin-3)/2;
                for i = 1:M
                    value = varargin{1*i};
                    str = varargin{2*i};
                    update_packet(obj, value, str);
                end
            else
                update_packet(obj, value, str);
            end
        end
    end
    methods (Access = private)
        function update_packet(obj, value, str)
            if any( contains(fieldnames(obj.ParaRecord),str) )
                obj.ParaRecord.(str) = value;
            else
                fprintf('Unknown keywords!\n');
            end
        end
    end
end