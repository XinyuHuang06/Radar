function varargout = general_normalize(varargin)
    input_num = varargin{1};
    input_num = abs(input_num);
    norm_value = norm(input_num,"fro");
    varargout{1} = input_num/norm_value;

end