% 参数设置
clear,clc
d = 0.003;       % 母泡直径 3mm (单位: m)
epsilon = 1.0;   % 能量耗散率 (m²/s³)
rho_c = 1000;    % 水密度 (kg/m³)
sigma = 0.0725;  % 表面张力 (N/m)
nu = 1e-6;       % 水运动粘度 (m²/s)

% 生成随机数
fv_sample = sampleDaughterSize(d, epsilon, rho_c, sigma, nu);

% 批量生成100个样本
n_samples = 1000;
samples = zeros(n_samples, 1);
tic
parfor i = 1:n_samples
    samples(i) = sampleDaughterSize(d, epsilon, rho_c, sigma, nu, ...
                                   'nPoints', 31); % 加速计算
end
toc
% 绘制直方图
histogram(samples,10, 'Normalization', 'pdf');
xlabel('Breakup fraction f_v');
ylabel('Probability density');
title('Daughter Size Distribution (d=3mm, \epsilon=1.0 m^2/s^3)');

function fv_sample = sampleDaughterSize(d, epsilon, rho_c, sigma, nu, varargin)
    % 输入参数:
    %   d: 母气泡直径 (m)
    %   epsilon: 能量耗散率 (m²/s³)
    %   rho_c: 连续相密度 (kg/m³)
    %   sigma: 表面张力 (N/m)
    %   nu: 连续相运动粘度 (m²/s)
    % 可选参数:
    %   'alpha_d': 分散相体积分数 (默认0.1)
    %   'delta': 防止奇异的参数 (默认0.01)
    %   'nPoints': f_v离散点数 (默认51)
    %   'lambda_ratio': λ_min/λ_η 比例 (默认31.4)
    %
    % 输出:
    %   fv_sample: 满足女儿尺寸分布的随机数 (f_v ∈ [0,1])
    
    % 解析可选参数
    p = inputParser;
    addParameter(p, 'alpha_d', 0.01, @isnumeric);
    addParameter(p, 'delta', 0.01, @isnumeric);
    addParameter(p, 'nPoints', 51, @isnumeric);
    addParameter(p, 'lambda_ratio', 31.4, @isnumeric);
    parse(p, varargin{:});
    
    alpha_d = p.Results.alpha_d;
    delta = p.Results.delta;
    nPoints = p.Results.nPoints;
    lambda_ratio = p.Results.lambda_ratio;
    
    % 计算Kolmogorov尺度和最小涡尺寸
    lambda_eta = (nu^3 / epsilon)^(1/4);    % Kolmogorov长度尺度
    lambda_min = lambda_ratio * lambda_eta;  % 最小有效涡尺寸
    lambda_max = d;                         % 最大涡尺寸（≤d）
    
    % 处理无效λ范围
    if lambda_min >= lambda_max
        fv_sample = 0.5;  % 返回中值作为默认
        return;
    end
    
    % 离散化f_v (0到0.5)
    fv_vec = linspace(1e-6, 0.5, nPoints)';  % 避免f_v=0
    
    % 计算每个f_v对应的积分值
    integral_vals = zeros(size(fv_vec));
    for i = 1:length(fv_vec)
        integral_vals(i) = computeIntegralForB(fv_vec(i), d, lambda_min, lambda_max, ...
                                             rho_c, epsilon, sigma, delta);
    end
    
    % 构建完整[0,1]区间分布（对称性）
    if mod(nPoints, 2) == 1
        mid_idx = (nPoints + 1)/2;
        fv_full = [fv_vec; 1 - fv_vec(1:mid_idx-1)];
        prob_full = [integral_vals; integral_vals(1:mid_idx-1)];
    else
        fv_full = [fv_vec; 0.5; 1 - fv_vec(1:end-1)];
        prob_mid = computeIntegralForB(0.5, d, lambda_min, lambda_max, ...
                                     rho_c, epsilon, sigma, delta);
        prob_full = [integral_vals; prob_mid; integral_vals(1:end-1)];
    end
    
    % 归一化概率
    total_prob = sum(prob_full);
    if total_prob <= 0
        fv_sample = 0.5;  % 返回中值作为默认
        return;
    end
    prob_full = prob_full / total_prob;
    
    % 生成随机数（离散逆变换采样）
    cdf = cumsum(prob_full);
    u = rand();
    idx = find(cdf >= u, 1);
    
    if isempty(idx)
        fv_sample = fv_full(end);
    elseif idx == 1
        fv_sample = fv_full(1);
    else
        % 线性插值
        fv_sample = fv_full(idx-1) + (fv_full(idx) - fv_full(idx-1)) * ...
                   (u - cdf(idx-1)) / (cdf(idx) - cdf(idx-1));
    end
end

% 辅助函数
function cf = c_f(fv)
    % 计算表面能增加系数
    cf = fv.^(2/3) + (1 - fv).^(2/3) - 1;
end

function fv_max = get_fv_max(cf_max)
    % 二分法求解f_v_max
    fv_low = 0;
    fv_high = 0.5;
    tol = 1e-6;
    
    if cf_max >= (2^(1/3) - 1)
        fv_max = 0.5;
        return;
    end
    
    while (fv_high - fv_low) > tol
        fv_mid = (fv_low + fv_high) / 2;
        cf_mid = c_f(fv_mid);
        
        if cf_mid > cf_max
            fv_high = fv_mid;
        else
            fv_low = fv_mid;
        end
    end
    fv_max = fv_mid;
end

function P_b = calc_P_b_given_e(fv, d, e_lambda, lambda, sigma, delta)
    % 计算P_b(f_v | d, e(λ), λ)
    % 使用元素级操作确保处理向量输入
    numerator = pi * lambda^3 * sigma;
    denominator = 6 * e_lambda * d;
    
    % 确保向量操作兼容
    fv_min = (numerator ./ denominator).^3;  % 使用元素级除法
    
    % 计算cf_max时也使用元素级操作
    cf_max = min(2^(1/3) - 1, e_lambda ./ (pi * d^2 * sigma));
    
    % 初始化输出
    P_b = zeros(size(e_lambda));
    
    for i = 1:numel(e_lambda)
        % 对每个e_lambda元素单独处理
        if cf_max(i) >= (2^(1/3) - 1)
            fv_max_i = 0.5;
        else
            fv_max_i = get_fv_max(cf_max(i));
        end
        
        if fv_min(i) >= fv_max_i || fv < fv_min(i) || fv > fv_max_i
            P_b(i) = 0;
        else
            P_b(i) = 1 / (fv_max_i - fv_min(i) + delta);
        end
    end
end

function P_b_lambda = calc_P_b_given_lambda(fv, d, lambda, rho_c, epsilon, sigma, delta)
    % 计算P_b(f_v | d, λ)
    mean_e = (pi/6) * rho_c * epsilon^(2/3) * lambda^(11/3);  % 平均涡动能
    
    if mean_e <= 0
        P_b_lambda = 0;
        return;
    end
    
    % 数值积分：e(λ)从0到20*mean_e
    e_upper = 20 * mean_e;
    
    % 定义被积函数 - 处理向量输入
    integrand = @(e) calc_P_b_given_e(fv, d, e, lambda, sigma, delta);
    
    P_b_lambda = integral(integrand, 0, e_upper, 'RelTol', 1e-3, 'AbsTol', 1e-6, 'ArrayValued', true);
end

function integral_val = computeIntegralForB(fv, d, lambda_min, lambda_max, ...
                                          rho_c, epsilon, sigma, delta)
    % 计算λ积分项 - 使用数组化计算确保输入为标量
    integrand = @(lambda) arrayfun(...
        @(lam) calc_P_b_given_lambda(fv, d, lam, rho_c, epsilon, sigma, delta) .* ...
               ((lam + d).^2 ./ lam.^(11/3)), ...
        lambda);
    
    integral_val = integral(integrand, lambda_min, lambda_max, ...
                          'RelTol', 1e-3, 'AbsTol', 1e-6, ...
                          'ArrayValued', true);
end