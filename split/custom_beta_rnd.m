function samples = custom_beta_rnd(a, b, n)
% 生成符合分布 y = (betapdf(x, a, b) + betapdf(1-x, a, b))/2 的随机数
% 输入:
%   a, b: Beta分布参数
%   n: 样本数量
% 输出:
%   samples: n×1向量，生成的随机数

% 生成标志位（决定每个样本来自哪个分布）
flag = rand(n, 1) < 0.5;

% 生成Beta(a, b)的样本
beta_samples = betarnd(a, b, n, 1);

% 将部分样本转换为1 - beta_samples（即Beta(b, a)的样本）
samples = beta_samples;
samples(~flag) = 1 - beta_samples(~flag); % 对非标志位样本取镜像
end