function [fit_coeff, integral_value, r_squared] = polynomial_fit_and_integral(x, y, degree, plot_results)
% 多项式拟合与积分计算
% 输入：
%   x: x坐标向量
%   y: y坐标向量
%   degree: 多项式阶数
%   plot_results: 是否绘制结果（1=是，0=否）
%
% 输出：
%   fit_coeff: 多项式系数（从高次到低次）
%   integral_value: 拟合曲线的积分值
%   r_squared: 拟合的R²决定系数

% 检查输入参数
if nargin < 4
    plot_results = 1; % 默认绘制结果
end

% 验证输入向量
if length(x) ~= length(y)
    error('错误：x和y向量长度必须相同！');
end
if degree >= length(x)
    warning('多项式阶数接近数据点数，可能导致过拟合！');
end

% 多项式拟合
[fit_coeff, ~] = polyfit(x, y, degree)

% 计算拟合值
y_fit = polyval(fit_coeff, x);

% 计算R²值
SS_res = sum((y - y_fit).^2);
SS_tot = sum((y - mean(y)).^2);
r_squared = 1 - (SS_res / SS_tot);

% 创建更密集的点用于平滑绘图
x_dense = linspace(min(x), max(x), 500);
y_dense = polyval(fit_coeff, x_dense);

% 计算拟合曲线的积分
integral_value = polyint(fit_coeff); % 求不定积分
integral_value = diff(polyval(integral_value, [min(x), max(x)])); % 计算定积分

% 显示结果
fprintf('=== 多项式拟合与积分结果 ===\n');
fprintf('多项式阶数: %d\n', degree);
fprintf('拟合多项式: ');
% poly_str = poly2str(fit_coeff, 'x');
% fprintf('%s\n', poly_str);
fprintf('R²决定系数: %.6f\n', r_squared);
fprintf('拟合曲线在[%.2f, %.2f]的积分值: %.6f\n', min(x), max(x), integral_value);

% 绘制结果
if plot_results
    figure;
    
    % 绘制原始数据和拟合曲线
    subplot(2, 1, 1);
    plot(x, y, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b'); % 原始数据点
    hold on;
    plot(x_dense, y_dense, 'r-', 'LineWidth', 2); % 拟合曲线
    hold off;
    title(sprintf('%d阶多项式拟合 (R² = %.4f)', degree, r_squared));
    xlabel('x');
    ylabel('y');
    legend('原始数据', '拟合曲线', 'Location', 'best');
    grid on;
    
    % 绘制积分区域
    subplot(2, 1, 2);
    plot(x_dense, y_dense, 'r-', 'LineWidth', 2); % 拟合曲线
    hold on;
    area(x_dense, y_dense, 'FaceColor', [0.8 0.9 1.0], 'EdgeColor', 'none'); % 填充积分区域
    plot(x, y, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b'); % 原始数据点
    hold off;
    title(sprintf('拟合曲线积分区域 (积分值 = %.4f)', integral_value));
    xlabel('x');
    ylabel('y');
    grid on;
    
    % 调整布局
    set(gcf, 'Position', [100, 100, 800, 600]);
end

% 返回结果
if nargout == 0
    clear fit_coeff integral_value r_squared
end
end