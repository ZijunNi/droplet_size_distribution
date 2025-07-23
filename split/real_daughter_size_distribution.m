clear,clc

d = 6*1.37e-2*1e-2;       % 气泡直径
re_tau = 180;
delta = 5e-3;
visc = 2.4e-6;
epsilon_est = re_tau^4*visc^3/delta^4;% 用近壁面visc*(du/dy)^2进行估计

epsilon = epsilon_est;% 24.76;  % 能量耗散率(Re_tau = 180 =>24.76 m²/s³)

tic
result = calculate_beta_fv(d, epsilon);
toc



% load("auto_saved_data_d_0.3_eps_24.76.mat");
% 创建示例数据
x =[result(1:end,1)];
y =[result(1:end,2)];

plot(x,y)
title(['Daughter Size Distribution ($d=',num2sci(d*1e+3,3),'$mm, ' ...
    '$\varepsilon=',num2sci(epsilon,3),' \rm{m^2/s^3}$)'],Interpreter='latex');
xlabel('$f_v$',Interpreter='latex')
ylabel('$\beta(f_v;d,\epsilon)$',Interpreter='latex')
%%
% 调用函数进行多项式拟合并计算积分
poly_coeff = polynomial_fit_and_integral(x, y, 22, 1);


% 生成10000个随机数
random_numbers = generate_poly_random(poly_coeff, 10000);

% 绘制直方图验证
histogram(random_numbers,20, 'Normalization', 'pdf');
hold on;

% 绘制理论密度曲线
x_plot = linspace(0, 0.5, 100);
pdf_norm = @(x) polyval(poly_coeff, x) / integral(@(x) polyval(poly_coeff, x), 0, 0.5);
plot(x_plot, pdf_norm(x_plot), 'r-', 'LineWidth', 2);
xlabel('x');
ylabel('Probability Density');
legend('Generated Samples', 'Theoretical PDF');

% 定义需要保存的变量名称（根据实际情况修改）
targetVars = {'d', 'epsilon', 'result'}; % 替换为你的变量名

fileName = ['auto_saved_data_d_', num2str(d*1000),'_eps_',num2str(epsilon) ,'.mat'];

save(fileName, targetVars{:});
% 显示保存信息
fprintf('变量已保存至文件: %s\n', fileName);
disp('保存变量列表:');
disp(targetVars);


%%
% 绘图
figure;
plot(result(:,1), result(:,2));
curve_integral(result(:,1),result(:,2))
xlabel('$f_v$',Interpreter='latex');
ylabel('$\beta(f_v,d)$',Interpreter='latex');
% xlim([0 0.1])
title(['Daughter Size Distribution ($d=',num2str(d*1e+3),'$mm, ' ...
    '$\epsilon=',num2str(epsilon),' \rm{m^2/s^3}$)'],Interpreter='latex');
grid on;



function beta_fv = calculate_beta_fv(d, epsilon, varargin)
% 计算气泡/液滴破碎的daughter size distribution beta(f_v, d)
% 参考文献：Wang et al. (2003). Chemical Engineering Science 58, 4629-4637
%
% 输入：
%   d        : 母气泡直径 (m)
%   epsilon  : 能量耗散率 (m^2/s^3)
% 可选参数：
%   sigma    : 表面张力 (N/m, 默认 0.0725 水-空气)
%   rho_c    : 连续相密度 (kg/m^3, 默认 1000 水)
%   nu       : 运动粘度 (m^2/s, 默认 1e-6 水)
%   delta    : 模型参数 (默认 0.01)
%   lambda_min_ratio : 最小涡尺寸比例 (默认 31.4)
%   n_lambda : 涡尺寸离散点数 (默认 100)
%   n_fv     : f_v离散点数 (默认 100)
%   n_e      : 涡动能离散点数 (默认 50)
%
% 输出：
%   beta_fv  : [f_v, beta] 矩阵，第一列为f_v，第二列为beta(f_v,d)

% ============== 参数设置 ==============
p = inputParser;
addParameter(p, 'sigma', 7.5e-3);    % 表面张力 (N/m)
addParameter(p, 'rho_c', 866);      % 密度 (kg/m^3)
addParameter(p, 'nu', 2.4e-6);         % 运动粘度 (m^2/s)
addParameter(p, 'delta', 0.05);      % 避免奇异的参数
addParameter(p, 'lambda_min_ratio', 1); % Kolmogorov尺度比例
addParameter(p, 'n_lambda', 400);    % λ离散点数
addParameter(p, 'n_fv', 400);        % f_v离散点数
addParameter(p, 'n_e', 400);          % e(λ)离散点数
parse(p, varargin{:});
param = p.Results;

% ============== 物理常数和离散化 ==============
% 验证Weber数是否大于临界值（We_crit≈1.2）
We = 2*param.rho_c * epsilon^(2/3) * d^(5/3) / param.sigma
if We < 1.2
    warning('Weber数%.2f < 1.2, 破碎可能不发生', We);
end

% 计算Kolmogorov长度尺度
eta = (param.nu^3 / epsilon)^(1/4); % Kolmogorov尺度 (m)
lambda_min = param.lambda_min_ratio * eta; % 最小涡尺寸 (m)
lambda_max = d; % 最大涡尺寸 (母气泡直径)
if(lambda_min>d)
    warning('无有效涡！lambda_min>d！')
end


% 离散化涡尺寸λ (对数分布)
lambda_eddy = logspace(log10(lambda_min), log10(lambda_max), param.n_lambda)';
d_lambda = diff(lambda_eddy);
d_lambda = [d_lambda; d_lambda(end)]; % 扩展最后一个差分

% 离散化破碎分数f_v (0.001 ~ 0.5)
f_v_half = linspace(0.001, 0.5, param.n_fv)';

% 初始化破碎核b(f_v|d) - 仅计算0-0.5区间
b_fv_half = zeros(size(f_v_half));

% ============== 预先计算c_f到f_v_max的映射表 ==============
% 创建插值表以提高计算速度
c_f_list = linspace(0, 2^(1/3)-1, 1000);
f_v_max_list = arrayfun(@(cf) solve_cf_equation(cf), c_f_list);

% ============== 主计算循环 ==============
for i = 1:length(lambda_eddy)
    lam = lambda_eddy(i);
    
    % 计算当前λ的平均涡动能 (公式13)
    u_lambda_mean = sqrt(2) * (epsilon * lam)^(1/3); % 平均涡速度 (m/s)
    e_lambda_mean = (pi/6) * lam^3 * param.rho_c * (u_lambda_mean^2)/2; % 平均动能 (J)
    
    % 离散化涡动能e(λ) (指数分布，从0到5倍平均动能)
    e_min = 1e-12 * e_lambda_mean; % 避免零
    e_max = 40 * e_lambda_mean; % 截断值（覆盖>99%的概率）
    e_eddy = linspace(e_min, e_max, param.n_e)';
    de = e_eddy(2) - e_eddy(1);
    
    % 计算涡动能分布 (公式12)
    P_e = (1/e_lambda_mean) * exp(-e_eddy/e_lambda_mean);
    
    % 计算碰撞频率 (公式5)
    omega_lambda = 0.923 * (1 - 0) * 1 * epsilon^(1/3) * (lam + d)^2 / lam^(11/3);
    
    for j = 1:length(f_v_half)
        f = f_v_half(j);
        
        % 计算当前f_v所需的最小能量 (公式6和9)
        % 毛细约束要求的最小动能
        d1_min = d * f^(1/3); % 子气泡直径
        e_capillary = pi * lam^3 * param.sigma / (6 * d1_min);
        
        % 表面能约束要求的最小动能
        c_f = f^(2/3) + (1-f)^(2/3) - 1;
        e_surface = c_f * pi * d^2 * param.sigma;
        
        % 综合约束：取两者最大值
        e_req = max(e_capillary, e_surface);
        
        % 找到满足e(λ) >= e_req的索引
        valid_e_idx = e_eddy >= e_req;
        
        if any(valid_e_idx)
            % 提取满足条件的涡动能
            e_valid = e_eddy(valid_e_idx);
            P_e_valid = P_e(valid_e_idx);
            
            % 初始化概率积分
            prob_integral = 0;
            
            for k = 1:length(e_valid)
                e_k = e_valid(k);
                
                % 计算当前e_k对应的边界条件 (公式7和10)
                % 最小破碎分数 (毛细压力约束)
                f_v_min_k = (pi * lam^3 * param.sigma / (6 * e_k * d))^3;
                
                % 最大破碎分数 (表面能约束)
                c_f_max_k = min(2^(1/3)-1, e_k / (pi * d^2 * param.sigma));
                f_v_max_k = interp1(c_f_list, f_v_max_list, c_f_max_k, 'linear', 'extrap');
                f_v_max_k = min(f_v_max_k, 0.5); % 确保≤0.5
                
                % 检查f_v是否在可行区间内
                if f >= f_v_min_k && f <= f_v_max_k
                    % 计算区间长度
                    L_interval = f_v_max_k - f_v_min_k;
                    
                    % 计算概率密度 (公式11)
                    prob_density = 1 / (L_interval + param.delta);
                    
                    % 累加到积分 (公式14)
                    prob_integral = prob_integral + prob_density * P_e_valid(k) * de;
                end
            end
            
            % 累加到破碎核 (公式1)
            b_fv_half(j) = b_fv_half(j) + prob_integral * omega_lambda * d_lambda(i);
        end
    end
end

% ============== 处理对称性和归一化 ==============
% 由于对称性，b(f_v|d) = b(1-f_v|d)，所以我们只需要计算0-0.5区间
% 但分母积分需要覆盖0-1整个区间

% 创建完整的f_v范围 (0-1)
f_v_full = [f_v_half; 1 - flip(f_v_half(1:end-1))]; % 避免重复0.5

% 创建完整的b(f_v|d)值 (对称复制)
b_fv_full = [b_fv_half; flip(b_fv_half(1:end-1))];

% 计算归一化因子 (公式17分母) - 整个[0,1]区间
integral_b_full = trapz(f_v_full, b_fv_full);

% 计算daughter size distribution (公式17)
% 注意：分子中的2倍因子是因为β(f_v,d)包含两个子气泡
beta_fv_half = 2 * b_fv_half / integral_b_full;

% 返回结果矩阵 (仅0-0.5区间，因为>0.5的部分是对称的)
beta_fv = [f_v_half, beta_fv_half];

% ============== 嵌套函数：求解f_v_max方程 ==============
function f_v_max = solve_cf_equation(c_f_max)
    % 求解方程: c_f = f^{2/3} + (1-f)^{2/3} - 1
    % 注意: f_v_max ∈ [0, 0.5]
    
    if c_f_max >= (2^(1/3)-1 - 1e-6)
        f_v_max = 0.5; % 最大可能值
    else
        % 数值求解方程 (使用fminbnd)
        fun = @(f) abs((f.^(2/3) + (1-f).^(2/3) - 1) - c_f_max);
        f_v_max = fminbnd(fun, 0, 0.5);
    end
end

end