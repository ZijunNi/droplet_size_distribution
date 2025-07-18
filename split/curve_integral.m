% 对一对向量进行梯形积分
function area = curve_integral(x, y)
% 检查输入
if nargin < 2
    error('需要输入x和y两个向量');
end
if length(x) ~= length(y)
    error('x和y长度不一致');
end

% 核心计算
area = trapz(x, y);
end