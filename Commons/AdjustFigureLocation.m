% Example: 
% :param : 
% :return : 
% detailed description: 
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 19/06/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
function AdjustFigureLocation()
    figs = findall(0, 'type', 'figure');
    screenSize = get(0, 'ScreenSize');
    screenWidth = screenSize(3);
    screenHeight = screenSize(4);
    numFigs = length(figs);
    numCols = ceil(sqrt(numFigs)); 
    numRows = ceil(numFigs / numCols);
    figWidth = screenWidth * 0.8 / numCols;
    figHeight = screenHeight * 0.8 / numRows;
    leftMargin = screenWidth * 0.1;
    topMargin = screenHeight * 0.1;
    for i = 1:numFigs
        row = ceil(i / numCols);
        col = mod(i-1, numCols) + 1;
        figLeft = leftMargin + (col-1) * figWidth;
        figBottom = screenHeight - topMargin - row * figHeight;
        set(figs(i), 'Position', [figLeft figBottom figWidth figHeight]);
    end
end