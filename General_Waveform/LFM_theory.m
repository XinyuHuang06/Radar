% Example: 
% :param : 
% :return : 
% detailed description: 
%------------------------------------------------------------------------------
% Created by: Xinyu Huang.
% On: 09/10/2024.
% Copyright (C) 2024 Xinyu Huang (learning_huang@163.com).
% All Rights Reserved.
% Unauthorized copying of this file, via any medium is strictly prohibited.
% Proprietary and confidential.
%------------------------------------------------------------------------------
 

fig1 = figure;
set(gcf, 'Units', 'centimeters', 'Position', [10, 5, 40, 12]); 
subplot(121);
hold on;
% 定义斜线的起始点和终止点
x1 = -5;
y1 = 2;
x2 = 5;
y2 = 6;
axis off
xlim([-10 10]);
ylim([-2 10]);
% 计算斜率
slope = (y2 - y1) / (x2 - x1);
% 绘制带箭头的坐标轴
quiver(-10, 0, 20, 0, 0, 'k', 'MaxHeadSize', 0.1,'LineWidth', 2); % X轴箭头
quiver(0, -5, 0, 15, 0, 'k', 'MaxHeadSize', 0.2,'LineWidth', 2); % Y轴箭头
% 绘制斜线
plot([x1, x2], [y1, y2], 'k-', 'LineWidth', 2);
% 虚线
plot([x1, x1], [y1, 0], 'k--', 'LineWidth', 2);
plot([x1, 0], [y1, y1], 'k--', 'LineWidth', 2);
plot([x2, x2], [y2, 0], 'k--', 'LineWidth', 2);
plot([x2, 0], [y2, y2], 'k--', 'LineWidth', 2);

text(x1,-1.2, '$-\frac{T}{2}$', 'Interpreter', 'latex', 'FontSize', 18,...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(x2,-1.2, '$\frac{T}{2}$', 'Interpreter', 'latex', 'FontSize', 18,...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(1.8, y1, '$f_0-\frac{B}{2}$', 'Interpreter', 'latex', 'FontSize', 18,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
text(-1.8, y2, '$f_0+\frac{B}{2}$', 'Interpreter', 'latex', 'FontSize', 18,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

text(10.5, 0, '$t$', 'Interpreter', 'latex', 'FontSize', 18,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
text(0, 10.5, '$f$', 'Interpreter', 'latex', 'FontSize', 18,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
% text(mid_x, mid_y+0.2, '$\mu$', 'Interpreter', 'latex', 'FontSize', 18);
% 标注斜率
mid_x = (x1 + x2) / 2;
mid_y = (y1 + y2) / 2;
text(mid_x, mid_y+0.2, '$\mu$', 'Interpreter', 'latex', 'FontSize', 18,...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
title('(a) 正斜率线性调频', 'Position', [0, -3, 0], 'HorizontalAlignment', 'center');

subplot(122);
hold on;
% 定义斜线的起始点和终止点
x1 = -5;
y1 = 6;
x2 = 5;
y2 = 2;
axis off
xlim([-10 10]);
ylim([-2 10]);
% 计算斜率
slope = (y2 - y1) / (x2 - x1);
% 绘制带箭头的坐标轴
quiver(-10, 0, 20, 0, 0, 'k', 'MaxHeadSize', 0.1,'LineWidth', 2); % X轴箭头
quiver(0, -5, 0, 15, 0, 'k', 'MaxHeadSize', 0.2,'LineWidth', 2); % Y轴箭头
% 绘制斜线
plot([x1, x2], [y1, y2], 'k-', 'LineWidth', 2);
% 虚线
plot([x1, x1], [y1, 0], 'k--', 'LineWidth', 2);
plot([x1, 0], [y1, y1], 'k--', 'LineWidth', 2);
plot([x2, x2], [y2, 0], 'k--', 'LineWidth', 2);
plot([x2, 0], [y2, y2], 'k--', 'LineWidth', 2);

text(x1,-1.2, '$-\frac{T}{2}$', 'Interpreter', 'latex', 'FontSize', 18,...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(x2,-1.2, '$\frac{T}{2}$', 'Interpreter', 'latex', 'FontSize', 18,...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(1.8, y1, '$f_0-\frac{B}{2}$', 'Interpreter', 'latex', 'FontSize', 18,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
text(-1.8, y2, '$f_0+\frac{B}{2}$', 'Interpreter', 'latex', 'FontSize', 18,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');

text(10.5, 0, '$t$', 'Interpreter', 'latex', 'FontSize', 18,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
text(0, 10.5, '$f$', 'Interpreter', 'latex', 'FontSize', 18,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
% text(mid_x, mid_y+0.2, '$\mu$', 'Interpreter', 'latex', 'FontSize', 18);
% 标注斜率
mid_x = (x1 + x2) / 2;
mid_y = (y1 + y2) / 2;
text(mid_x, mid_y+0.2, '$\mu$', 'Interpreter', 'latex', 'FontSize', 18,...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

title('(b) 负斜率线性调频', 'Position', [0, -3, 0], 'HorizontalAlignment', 'center');

exportgraphics(fig1, './output/LFM_theory.jpg', 'Resolution', 300);
% saveas(fig1, './output/LFM_theory.jpg')