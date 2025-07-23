function samples = custom_rand_ao_shape(n, a, h_low)
% 生成分段均匀的"凹"字形分布随机数
% 输入:
%   n     : 样本数量
%   a     : 左侧区间结束点 (0 < a < 0.5)
%   h_low : 中间区间(a,1-a)的概率密度
% 输出:
%   samples : 生成的随机数向量 (n×1)

% 参数验证
if a <= 0 || a >= 0.5
    error('a必须在(0,0.5)区间内');
end

% 计算右侧区间开始点(对称)
b = 1 - a;

% 计算两侧区域概率密度 (归一化条件)
h_high = (1 - (1 - 2*a)*h_low) / (2*a);

% 验证概率密度有效性
if h_high <= 0
    error('h_low值过大导致h_high<=0! 最大允许h_low=%.4f', 1/(1-2*a));
end

% 计算CDF关键点
Fa = a * h_high;         % CDF(a)
Fb = Fa + (b - a)*h_low; % CDF(b)

% 生成均匀分布随机数
u = rand(n, 1);

% 预分配输出数组
samples = zeros(n, 1);

% 分段处理 (向量化计算提高效率)
% 左侧区间 [0, a]
left_mask = (u <= Fa);
samples(left_mask) = u(left_mask) / h_high;

% 中间区间 (a, b]
mid_mask = (u > Fa) & (u <= Fb);
samples(mid_mask) = a + (u(mid_mask) - Fa) / h_low;

% 右侧区间 (b, 1]
right_mask = (u > Fb);
samples(right_mask) = b + (u(right_mask) - Fb) / h_high;
end