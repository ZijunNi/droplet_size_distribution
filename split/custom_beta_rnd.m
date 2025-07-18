function samples = custom_beta_rnd(a, b, n)
% 生成符合分布 y = (betapdf(x, a, b) + betapdf(1-x, a, b))/2 的随机数
% 输入:
%   a, b: Beta分布参数
%   n: beta分布占比，0.1左右较佳
% 输出:
%   samples: n×1向量，生成的随机数

flag = 1;
sigma = 0;% 均匀分布参数
ratio = n;
while(flag)%随机数在0-1之间
    u = rand;
    if u < ratio
        % 从正态分布 N(0.5, 1) 生成样本
        % samples = unifrnd(sigma, 1-sigma);
        samples = normrnd(0.5, 0.3);
    else
        % 从两个Beta分布的混合生成样本
        if rand < 0.5
            % Beta(alpha, beta)
            samples = betarnd(a, b);
        else
            % Beta(beta, alpha) 即 1 - X, X ~ Beta(alpha, beta)
            samples = 1 - betarnd(a, b);
        end
    end
    if(abs(samples-0.5)<0.5)
        flag = 0;
    end
end


% %%%以下为旧beta模型，仅有beta分布
%     % 生成标志位（决定每个样本来自哪个分布）
%     flag = rand(n, 1) < 0.5;
% 
%     % 生成Beta(a, b)的样本
%     beta_samples = betarnd(a, b, n, 1);
% 
%     % 将部分样本转换为1 - beta_samples（即Beta(b, a)的样本）
%     samples = beta_samples;
%     samples(~flag) = 1 - beta_samples(~flag); % 对非标志位样本取镜像
% %%%旧模型结束

end