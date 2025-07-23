function [abs_log_residuals, in_range] = calcLogResiduals(x_ref, y_ref, x_target, y_target)
% CALCLOGRESIDUALS 计算离散点与参考曲线的对数空间绝对残差
%   [abs_log_residuals, in_range] = calcLogResiduals(x_ref, y_ref, x_target, y_target)
%
% 输入参数:
%   x_ref, y_ref - 参考曲线的离散点坐标 (向量)
%   x_target, y_target - 目标离散点坐标 (向量)
%
% 输出参数:
%   abs_log_residuals - 对数空间的绝对残差 (|log10(y_target) - log10(y_interp)|)
%   in_range - 逻辑数组，指示哪些目标点在参考曲线的x范围内
%
% 说明:
%   该函数在对数空间计算目标点与参考曲线之间的垂直距离残差，
%   并返回不带符号的绝对残差值

% ===== 数据预处理 =====
% 确保输入为列向量
x_ref = x_ref(:);
y_ref = y_ref(:);
x_target = x_target(:);
y_target = y_target(:);

% 1. 确保所有Y值为正（对数计算的前提）
y_ref(y_ref <= 0) = eps;
y_target(y_target <= 0) = eps;

% 2. 移除NaN或Inf值
valid_ref = isfinite(x_ref) & isfinite(y_ref);
x_ref = x_ref(valid_ref);
y_ref = y_ref(valid_ref);

valid_target = isfinite(x_target) & isfinite(y_target);
x_target = x_target(valid_target);
y_target = y_target(valid_target);

% 3. 确保参考曲线是单调递增的（对于插值很重要）
[x_ref, idx] = sort(x_ref);
y_ref = y_ref(idx);

% 在对数空间计算残差
log_y_ref = log10(y_ref);
log_y_target = log10(y_target);

% ===== 安全插值处理 =====
% 确定内插范围
xmin = min(x_ref);
xmax = max(x_ref);

% 只对范围内的点进行插值
in_range = (x_target >= xmin) & (x_target <= xmax);
x_target_inrange = x_target(in_range);

% 初始化插值结果（全NaN）
log_y_interp = NaN(size(x_target));

% 对范围内的点进行插值
if ~isempty(x_target_inrange)
    % 使用pchip插值方法
    log_y_interp(in_range) = interp1(x_ref, log_y_ref, x_target_inrange, 'pchip');
end

% 计算对数空间的残差（带符号）
log_residuals = log_y_target - log_y_interp;

% 计算绝对残差（不带符号）
abs_log_residuals = abs(log_residuals);
end