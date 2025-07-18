clear,clc
% 使用示例：将物理参数转换为无量纲参数
%%% Wang Model Ref.
% d = 3e-3;       % 5mm 气泡
% epsilon = 1;  % 能量耗散率 (m²/s³)
% rho_c = 1000;   % 水密度 (kg/m³)
% sigma = 0.07;  % 表面张力 (N/m)
% nu = 1e-6;      % 运动粘度 (m²/s)
%%% 以上给出We=0.8915, nu_star = 0.0023

%%% Exp. Ref
d = 3e-4;       % 5mm 气泡
epsilon = 10000;  % 能量耗散率 (m²/s³)
rho_c = 1000;   % 水密度 (kg/m³)
sigma = 0.0057;  % 表面张力 (N/m)
nu = 2.4e-6;      % 运动粘度 (m²/s)
tic
[We, nu_star] = to_dimensionless(d, epsilon, rho_c, sigma, nu);

%%

% 计算无量纲分布
result = calculate_beta_fv_dimensionless(We, nu_star);
toc
% 绘图
figure;
plot(result(:,1), result(:,2));
xlabel('破碎分数 f_v');
ylabel('概率密度 \beta(f_v,d)');
% xlim([0 0.1])
title(['Daughter Size Distribution (d=',num2str(d*1e+3),'mm, \epsilon=',num2str(epsilon),' m^2/s^3)']);
grid on;



function [We, nu_star] = to_dimensionless(d, epsilon, rho_c, sigma, nu)
% 将物理参数转换为无量纲参数
% 输入：
%   d       : 气泡直径 (m)
%   epsilon : 能量耗散率 (m²/s³)
%   rho_c   : 连续相密度 (kg/m³)
%   sigma   : 表面张力 (N/m)
%   nu      : 运动粘度 (m²/s)
% 输出：
%   We      : 韦伯数
%   nu_star : 无量纲粘度

We = (rho_c * epsilon^(2/3) * d^(5/3)) / sigma;
nu_star = nu / (epsilon^(1/3) * d^(4/3));

end

function beta_fv = calculate_beta_fv_dimensionless(We, nu_star, varargin)
% 无量纲形式计算daughter size distribution
% 输入：
%   We       : 韦伯数
%   nu_star  : 无量纲粘度
% 可选参数：
%   delta    : 模型参数 (默认 0.01)
%   lambda_min_ratio : 最小涡尺寸比例 (默认 31.4)
%   n_lambda : 涡尺寸离散点数 (默认 100)
%   n_fv     : f_v离散点数 (默认 100)
%   n_e      : 涡动能离散点数 (默认 50)
%
% 输出：
%   beta_fv  : [f_v, beta] 矩阵

% ============== 参数设置 ==============
p = inputParser;
addParameter(p, 'delta', 0.01);      % 避免奇异的参数
addParameter(p, 'lambda_min_ratio', 31.4); % Kolmogorov尺度比例
addParameter(p, 'n_lambda', 100);    % λ*离散点数
addParameter(p, 'n_fv', 100);        % f_v离散点数
addParameter(p, 'n_e', 50);          % e*离散点数
parse(p, varargin{:});
param = p.Results;

% ============== 计算无量纲最小涡尺寸 ==============
eta_star = nu_star^(3/4); % 无量纲Kolmogorov尺度
lambda_min_star = param.lambda_min_ratio * eta_star; % 最小涡尺寸
lambda_max_star = 1; % 最大涡尺寸 (母气泡直径)

% 离散化无量纲涡尺寸λ* (对数分布)
lambda_star = logspace(log10(lambda_min_star), log10(lambda_max_star), param.n_lambda)';
d_lambda_star = diff(lambda_star);
d_lambda_star = [d_lambda_star; d_lambda_star(end)]; % 扩展最后一个差分

% 离散化破碎分数f_v (0.001 ~ 0.5)
f_v_half = linspace(0.001, 0.5, param.n_fv)';

% 初始化破碎核b*(f_v)
b_fv_half = zeros(size(f_v_half));

% ============== 预先计算c_f到f_v_max的映射表 ==============
c_f_list = linspace(0, 2^(1/3)-1, 1000);
f_v_max_list = arrayfun(@(cf) solve_cf_equation(cf), c_f_list);

% ============== 主计算循环 ==============
for i = 1:length(lambda_star)
    lam_star = lambda_star(i);
    
    % 计算当前λ*的平均无量纲涡动能
    e_mean_star = (pi/6) * lam_star^(11/3);
    
    % 离散化无量纲涡动能e*
    e_min_star = 1e-12 * e_mean_star;
    e_max_star = 5 * e_mean_star;
    e_eddy_star = linspace(e_min_star, e_max_star, param.n_e)';
    de_star = e_eddy_star(2) - e_eddy_star(1);
    
    % 计算涡动能分布
    P_e = (1/e_mean_star) * exp(-e_eddy_star/e_mean_star);
    
    % 计算无量纲碰撞频率
    omega_lambda_star = 0.923 * (lam_star + 1)^2 / lam_star^(11/3);
    
    for j = 1:length(f_v_half)
        f = f_v_half(j);
        
        % 计算当前f_v的c_f值
        c_f = f^(2/3) + (1-f)^(2/3) - 1;
        
        % 初始化概率积分
        prob_integral = 0;
        
        for k = 1:length(e_eddy_star)
            e_star = e_eddy_star(k);
            
            % 计算毛细约束 (无量纲形式)
            f_v_min_k = (pi * lam_star^3 / (6 * e_star * We))^3;
            
            % 计算表面能约束 (无量纲形式)
            c_f_max_k = min(2^(1/3)-1, e_star * We / pi);
            f_v_max_k = interp1(c_f_list, f_v_max_list, c_f_max_k, 'linear', 'extrap');
            f_v_max_k = min(f_v_max_k, 0.5); % 确保≤0.5
            
            % 检查f_v是否在可行区间内
            if f >= f_v_min_k && f <= f_v_max_k
                % 计算区间长度
                L_interval = f_v_max_k - f_v_min_k;
                
                % 计算概率密度
                prob_density = 1 / (L_interval + param.delta);
                
                % 累加到积分
                prob_integral = prob_integral + prob_density * P_e(k) * de_star;
            end
        end
        
        % 累加到破碎核
        b_fv_half(j) = b_fv_half(j) + prob_integral * omega_lambda_star * d_lambda_star(i);
    end
end

% ============== 处理对称性和归一化 ==============
% 创建完整的f_v范围 (0-1)
f_v_full = [f_v_half; 1 - flip(f_v_half(1:end-1))];

% 创建完整的b(f_v)值 (对称复制)
b_fv_full = [b_fv_half; flip(b_fv_half(1:end-1))];

% 计算归一化因子
integral_b_full = trapz(f_v_full, b_fv_full);

% 计算daughter size distribution
beta_fv_half = 2 * b_fv_half / integral_b_full;

% 返回结果矩阵
beta_fv = [f_v_half, beta_fv_half];

% ============== 嵌套函数：求解f_v_max方程 ==============
function f_v_max = solve_cf_equation(c_f_max)
    if c_f_max >= (2^(1/3)-1 - 1e-6)
        f_v_max = 0.5;
    else
        fun = @(f) abs((f.^(2/3) + (1-f).^(2/3) - 1) - c_f_max);
        f_v_max = fminbnd(fun, 0, 0.5);
    end
end

end