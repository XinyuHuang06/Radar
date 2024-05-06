classdef RadarSignal < handle
    properties(Access = private)
        % Private Member
        s_N;
        s_fs;
        
    end
    properties(Access = public)
        % Public Member
        MatFilOut;
        PSL;
        ISL;
        PSLR;
        ISLR;
        PAPR;
        AmbiFun;
        CycSpec;
        s_name = 'none';
        signal;
    end
    methods
        % Construct Function
        function obj = RadarSignal(signal,name)
            obj.signal = signal;
            obj.s_N = length(signal);
            obj.s_name = name;
            obj.PAPR = Evaluate_PAPR(obj.signal);
            [obj.MatFilOut ,obj.PSL, obj.ISL, obj.PSLR, obj.ISLR] = Analysis_Sidelobe(obj.signal,obj.signal);
        end
    end
end