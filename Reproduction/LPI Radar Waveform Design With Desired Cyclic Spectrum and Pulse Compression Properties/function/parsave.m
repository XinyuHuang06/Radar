function parsave(fname, varargin)
    s = struct();
    for i = 1:numel(varargin)
        s.(inputname(1+i)) = varargin{i};
    end
    save(fname, '-struct', 's');
end

