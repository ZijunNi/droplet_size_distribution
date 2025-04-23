% A.m 文件内容
% Call ALL Function
function F = distribution_functions
F.symmetric_gaussian = @symmetric_gaussian;
F.symmetric_gamma = @symmetric_gamma;
F.symmetric_beta = @symmetric_beta;
F.skewed_double_beta = @skewed_double_beta;
end


function beta = symmetric_gaussian(f_v, mu, sigma, C)
    % 对称双高斯分布模型
    % 输入参数:
    %   f_v   : 体积分数 (0~1)
    %   mu    : 高斯峰中心位置 (需满足 mu < 0.5)
    %   sigma : 高斯峰标准差
    %   C     : 幅度缩放系数
    % 输出:
    %   beta  : 分布函数值
    
    % 计算双高斯分布
    gauss_left  = exp(-(f_v - mu).^2 / (2*sigma^2));
    gauss_right = exp(-(f_v - (1-mu)).^2 / (2*sigma^2));
    beta = C * (gauss_left + gauss_right);% .* f_v .* (1 - f_v);
    
    % 端点强制归零
    beta(f_v == 0 | f_v == 1) = 0;
end

function beta = symmetric_gamma(f_v, k, theta, C)
    % 对称修正伽马分布模型
    % 输入参数:
    %   f_v   : 体积分数 (0~1)
    %   k     : 伽马分布形状参数 (k > 1 时在 f_v=0 附近形成峰)
    %   theta : 伽马分布尺度参数
    %   C     : 幅度缩放系数
    % 输出:
    %   beta  : 分布函数值
    
    % 计算伽马分布及其镜像
    gamma_left  = (f_v.^(k-1) .* exp(-f_v/theta)) / (gamma(k)*theta^k);
    gamma_right = ((1 - f_v).^(k-1) .* exp(-(1 - f_v)/theta)) / (gamma(k)*theta^k);
    beta = C * (gamma_left + gamma_right) ;%.* f_v .* (1 - f_v);
    
    % 处理数值稳定性 (避免除以零)
    beta(isnan(beta)) = 0;
    beta(f_v == 0 | f_v == 1) = 0;
end

function beta = symmetric_beta(f_v, alpha,beta, C)
    % 对称修正贝塔分布模型
    % 输入参数:
    %   f_v   : 体积分数 (0~1)
    %   alpha : 贝塔分布形状参数 (alpha < 1 时增强双峰特性)
    %   C     : 幅度缩放系数
    % 输出:
    %   beta  : 分布函数值
    
    % 计算对称贝塔分布
    beta_pdf = betapdf(f_v, alpha, beta);
    beta = C * beta_pdf;% .* f_v .* (1 - f_v);
    
    % 端点强制归零
    % beta(f_v == 0 | f_v == 1) = 0;
end

function beta = skewed_double_beta(f_v, alpha, gamma, C)
    % 对称双峰Beta分布模型，支持峰偏移和中间非零
    % 输入参数:
    %   f_v   : 体积分数 (0~1)
    %   alpha : 控制峰尖锐度 (alpha < 1时峰更尖锐)
    %   gamma : 控制峰偏移程度 (0 < gamma < 1, 越大峰越靠近端点)
    %   C     : 整体幅值缩放
    % 输出:
    %   beta  : 分布函数值
    
    % 计算峰偏移位置参数
    mu = 0.5 * (1 - gamma);  % 峰中心位置 (gamma=0时mu=0.5，gamma=1时mu=0)
    
    % 左峰和右峰的Beta分布参数
    beta_left  = betapdf(f_v, alpha, (1 - mu) * alpha / mu);
    beta_right = betapdf(f_v, (1 - mu) * alpha / mu, alpha);
    
    % 组合双峰并乘窗口函数
    beta = C * (beta_left + beta_right) .* (f_v .* (1 - f_v)).^(1-alpha);
    
    % 端点归零
    beta(f_v == 0 | f_v == 1) = 0;
end