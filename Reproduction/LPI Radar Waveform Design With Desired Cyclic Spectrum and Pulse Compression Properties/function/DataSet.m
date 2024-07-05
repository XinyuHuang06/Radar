
classdef DataSet < handle
    properties (Access = public)
        packets;
    end
    methods
        function obj = DataSet(xr, br, cr, lambda_0, lambda_1, rho_0, rho_1, r, h, vartheta, N)
            obj.packets.xr = xr;
            obj.packets.br = br;
            obj.packets.cr = cr;
            obj.packets.rho_0 = rho_0;
            obj.packets.rho_1 = rho_1;
            obj.packets.lambda_0 = lambda_0;
            obj.packets.lambda_1 = lambda_1;
            obj.packets.r = r;
            obj.packets.h = h;
            obj.packets.vartheta = vartheta;
            obj.packets.N = N;
        end
    end
    methods (Access = public)
        function update(obj, value, str, varargin)
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
            if any( contains(fieldnames(obj.packets),str) )
                obj.packets.(str) = value;
            else
                fprintf('Unknown keywords!\n');
            end
        end
    end
end



