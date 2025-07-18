function slope = logLinearFit(x, y)
% LOGLINEARFIT 对y取对数后与x进行线性拟合并返回斜率
%   输入：
%       x: 横坐标数据向量
%       y: 纵坐标数据向量（需全为正数）
%   输出：
%       slope: 拟合直线的斜率

% 检查输入有效性
if nargin < 2
    error('需要两个输入参数：x和y');
end
if ~isvector(x) || ~isvector(y)
    error('x和y必须是向量');
end
if length(x) ~= length(y)
    error('x和y的长度必须相同');
end
if any(y <= 0)
    error('y必须全为正数');
end

% 对y取自然对数
logy = log(y);

% 线性拟合（使用多项式拟合的一阶形式）
p = polyfit(x, logy, 1);

% 提取斜率（一阶项系数）
slope = p(1);
end