function random_numbers = generate_poly_random(poly_coeff, num_samples)
    % poly_coeff: 多项式系数向量（从高阶到低阶），例如 [a_n, a_{n-1}, ..., a_0]
    % num_samples: 生成的随机数数量
    
    % 定义积分区间
    a = 0;
    b = 0.5;
    
    % 创建多项式函数句柄
    poly_fun = @(x) polyval(poly_coeff, x);
    
    % 计算归一化常数（多项式在[a,b]上的积分）
    C = integral(poly_fun, a, b);
    
    % 归一化概率密度函数
    pdf_norm = @(x) poly_fun(x) / C;
    
    % 计算累积分布函数（CDF）的逆函数（数值方法）
    N = 1000; % 离散化点数
    x_vals = linspace(a, b, N);
    pdf_vals = pdf_norm(x_vals);
    cdf_vals = cumtrapz(x_vals, pdf_vals); % 数值积分计算CDF
    
    % 确保CDF单调递增（避免数值误差）
    [cdf_vals, unique_idx] = unique(cdf_vals);
    x_vals = x_vals(unique_idx);
    
    % 创建CDF逆函数的插值函数
    inv_cdf_fun = @(u) interp1(cdf_vals, x_vals, u, 'pchip');
    
    % 生成均匀分布的随机数并转换
    u = rand(num_samples, 1); % [0,1]均匀分布
    random_numbers = inv_cdf_fun(u);
end