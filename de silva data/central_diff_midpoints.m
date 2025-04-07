function [dy_dx, x_mid] = central_diff_midpoints(y, x)
    % 使用中心差分法计算中点处的导数值
    %
    % 输入参数:
    %   y : 函数值数组 (1×N 或 N×1)
    %   x : 对应的空间坐标 (1×N 或 N×1)
    %
    % 输出参数:
    %   dy_dx : 中点处的导数数组 (长度比y少1)
    %   x_mid : 中点坐标数组 (长度比x少1)

    % 确保输入是列向量
    y = y(:);
    x = x(:);
    
    % 检查输入长度是否匹配
    if length(y) ~= length(x)
        error('输入数组y和x的长度必须相同');
    end
    
    % 计算中点坐标
    x_mid = (x(1:end-1) + x(2:end)) / 2;
    
    % 计算中点处的导数
    dy_dx = diff(y) ./ diff(x);
end