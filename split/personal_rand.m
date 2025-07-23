function rand_num = personal_rand(i)
%自定义分布随机数生成

% 步骤1：定义PDF的离散数据点
x = linspace(0, 1, 1000); % 覆盖足够范围
pdf_values = abs(0.5-x); % 目标PDF

% 步骤2：归一化PDF（此处已归一化，但建议保留此步骤）
dx = x(2) - x(1);
pdf_normalized = pdf_values / (sum(pdf_values) * dx); % 确保积分归一化

% 步骤3：计算CDF
cdf = cumtrapz(x, pdf_normalized);
cdf = cdf / cdf(end); % 确保CDF末端为1

% 步骤4：生成均匀分布的随机数
N = 1e5; % 样本数量
u = rand(N, 1);

% 步骤5：插值获得对应x值（使用pchip保持单调性）
x_samples = interp1(cdf, x, u, 'pchip');

% 步骤6：处理可能的超出范围的值
x_samples = max(x_samples, x(1)); % 确保不低于最小x
x_samples = min(x_samples, x(end)); % 确保不高于最大x

end