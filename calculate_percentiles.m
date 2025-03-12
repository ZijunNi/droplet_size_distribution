function [q5, q95] = calculate_percentiles(vector,percent)
    % 计算向量的前 5% 分位和 95% 分位的值，并向中间舍入
    %
    % 输入：
    %   vector - 输入的长向量
    %
    % 输出：
    %   q5  - 前 5% 分位的值（向中间舍入）
    %   q95 - 前 95% 分位的值（向中间舍入）
    
    % 检查输入是否为空
    if isempty(vector)
        error('输入向量不能为空。');
    end
    
    % 计算分位值
    q5 = quantile(vector, percent);  
    q95 = quantile(vector, 1-percent); 
    
  
end