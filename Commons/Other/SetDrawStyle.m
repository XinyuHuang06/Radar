function SetDrawStyle(fig)
% Example:
% :param fig: figure
% :return :
%---------------
% Created by: Xinyu Huang.
% On: 09/06/2023.
% Copyright (C) 2023 Xinyu.Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%---------------
    if nargin < 1
        fig = gcf; 
    end
    % Set the figure position
    screen_size = get(0,'ScreenSize'); % Getting the current screen size
    height = 600; % Set the height, width , the left and botton positions for the figure
    width = 800;
    left = screen_size(3)/2 - width/2;
    bottom = screen_size(4)/2 - height/2;
    set(fig, 'Position',[left,bottom,width,height]);
    % Set the font size
    TickLabelFontSize = 12; % Set the tick label font size
    AxisLabelFontSize = 14; % Set the axis label font size
    LegendFontSize = 8; % Set the legend font size
    TitleFontSize = 14; % Set the title font size
    FontWeight = 'bold'; % Set the font weight
    % The list of handles
    h_list = get(fig, 'children');
    % Loop over each handle
    for iHandle = 1:length(h_list)
        % If an axis, apply settings
        if(strcmp(get(h_list(iHandle), 'type'), 'axes'))
            h_child = get(h_list(iHandle), 'children');
            for iChild = 1:length(h_child)
                if strcmp(get(h_child(iChild), 'type'), 'line')
                    set(h_child(iChild), 'linewidth', 2)
                end
            end
            set(h_list(iHandle), 'fontweight', FontWeight)
            set(h_list(iHandle), 'fontsize', TickLabelFontSize)
            set(h_list(iHandle), 'fontname', 'Times New Roman')
            %  title settings
            h = get(h_list(iHandle), 'title');
            set(h, 'fontweight', FontWeight)
            set(h, 'fontsize', TitleFontSize)
            
            try % xlabel settings
                h = get(h_list(iHandle), 'xlabel');
                set(h, 'fontweight', FontWeight)
                set(h, 'fontsize', AxisLabelFontSize)
            catch
            end

            try % ylabel settings
                h = get(h_list(iHandle), 'ylabel');
                set(h, 'fontweight', FontWeight)
                set(h, 'fontsize', AxisLabelFontSize)
            catch
            end
            
            try % zlabel settings
                h = get(h_list(iHandle), 'zlabel');
                set(h, 'fontweight', FontWeight)
                set(h, 'fontsize', AxisLabelFontSize)
            catch
            end
            
            try % legend settings
                h = getget(h_list(iHandle), 'Legend');
                set(h, 'fontsize', LegendFontSize);
            catch  
            end
        end
    end
    ax = findobj(fig,'Type','axes'); 
    if ~isempty(ax)
        for iter_i = 1:length(ax)
            axes = ax(iter_i);
            set(axes, 'LooseInset', [0,0,0,0]);
        end
    end
end