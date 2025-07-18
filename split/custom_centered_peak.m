function samples = custom_centered_peak(b,num_samples)
    % 使用变换的拉普拉斯分布（更接近背靠背指数）
    % b = 0.7;  % 尺度参数（控制集中程度）
    u = rand(num_samples, 1) - 0.5;
    x = -b * sign(u) .* log(1 - 2*abs(u));  % 生成拉普拉斯分布
    samples = 1./(1 + exp(-x));             % 使用logistic函数压缩到(0,1)
end