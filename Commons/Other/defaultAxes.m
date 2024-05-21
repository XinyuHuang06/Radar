function defaultAxes
    ax = gca; hold on; box on
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.LineWidth = 1.2;
    ax.GridLineStyle = ':';
    ax.FontName = 'Cambria';
    ax.FontSize = 12;
    ax.GridAlpha = .5;
end