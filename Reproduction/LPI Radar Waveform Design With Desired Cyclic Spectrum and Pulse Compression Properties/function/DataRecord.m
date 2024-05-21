classdef DataRecord < handle
    properties (Access = public)
        DataBag;
    end
    methods
        function obj = DataRecord(Max_ItersNum)
            obj.DataBag.SUM = zeros(Max_ItersNum+1,1);
            obj.DataBag.TargetFun = zeros(Max_ItersNum+1,1);
            obj.DataBag.TargetLagrange_1 = zeros(Max_ItersNum+1,1);
            obj.DataBag.TargetLagrange_2 = zeros(Max_ItersNum+1,1);
            obj.DataBag.TargetPenalty_1 = zeros(Max_ItersNum+1,1);
            obj.DataBag.TargetPenalty_2 = zeros(Max_ItersNum+1,1);
        end
    end

    methods (Access = public)
        function UpdateDataSet(obj, DataBag, num)
            obj.DataBag.SUM(num) = DataBag.SUM;
            obj.DataBag.TargetFun(num) = DataBag.TargetFun;
            obj.DataBag.TargetLagrange_1(num) = DataBag.TargetLagrange_1;
            obj.DataBag.TargetLagrange_2(num) = DataBag.TargetLagrange_2;
            obj.DataBag.TargetPenalty_1(num)= DataBag.TargetPenalty_1;
            obj.DataBag.TargetPenalty_2(num) = DataBag.TargetPenalty_2;
        end
    end
end