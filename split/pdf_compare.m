
% Re = [100 180 200 400 800 1600 3200];
% uinf = [16.512925 18.2710 18.245793  19.978661 21.711529 23.444397 25.177265];
% semilogx(Re,uinf,'x');
% hold on
% semilogx(Re,2.5*(log(Re)-log(180))+uinf(2));
%
clear,clc;
% load("converted_de_silva_data.mat")
load("converted_de_silva_data.mat")
loc = 2;
figure;imagesc(xpos_delta,ypos_delta,squeeze(U(:,:,loc))'); %daspect([1 1 1]); %  wall parallel plane 
% U(:,:,1) = 0;
% save("converted_de_silva_data.mat","U","V","W","xpos_delta","ypos_delta",'zpos_delta')
% ylim([0 32])
a = squeeze(U(:,:,loc))';
mean(a(:))

%%

alpha = 2;      % 替换为所需的alpha值
beta = 15;     % 替换为所需的beta值
num_samples = 10000; % 样本数量
sigma = 0;
ratio = 0;

samples_beta = zeros(num_samples, 1);
samples_uniform = zeros(num_samples, 1);

i = 1;
while(i <= num_samples)
    samples_beta(i) = custom_beta_rnd(alpha, beta, ratio);
    samples_uniform(i) = custom_rand_ao_shape(1, 0.15, 0.2);
    i = i + 1;
end

f_v = linspace(0, 1, 1000);  
pdf = (1-ratio) * (betapdf(f_v, alpha, beta) + betapdf(1 - f_v, alpha, beta))/2 ...
      + ratio * unifpdf(f_v, sigma, 1-sigma);
pdf = 2*pdf;% 满足[0 0.5]半区间积分为1

figure;
% histogram(samples_beta,6, 'Normalization', 'pdf', 'DisplayName', '混合beta分布');
hold on;
load("no_important.mat")
% plot(f_BV, eta_v, 'b-', 'LineWidth', 2, 'DisplayName', 'Luo Model');
% histogram(samples_uniform,6, 'Normalization', 'pdf','DisplayName','混合均匀分布');
plot(f_v, pdf, 'r', 'LineWidth', 2, 'DisplayName', '理论PDF');
xlabel('f_v');
ylabel('PDF');
% set(gca,'YScale','log')
legend();
title('混合Beta分布样本 vs 理论PDF');
hold off;


%%
% 步骤1：定义PDF的离散数据点
x = linspace(0, 1, 1000); % 覆盖足够范围

k = 4;%k>4
a = 1/sqrt(k);
% if (x>0.5-a&&x<0.5+a)
%     pdf_values = abs(k*(0.5-x)); % 目标PDF
% else
%     pdf_values = 0;
% end
pdf = @(x) ((x > 0.5 - a) & (x < 0.5 + a)) .* k .* abs(0.5 - x);
pdf_values = pdf(x);
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

% % 步骤6：处理可能的超出范围的值
% x_samples = max(x_samples, x(1)); % 确保不低于最小x
% x_samples = min(x_samples, x(end)); % 确保不高于最大x

% 验证结果：绘制直方图与目标PDF对比
figure;
histogram(x_samples, 'Normalization', 'pdf', 'BinWidth', 0.2);
hold on;
plot(x, pdf_normalized, 'r', 'LineWidth', 2);
xlabel('x');
ylabel('PDF');
legend('Generated Samples', 'Target PDF');
title('Generated Samples vs Target PDF');
%%

% v1 = (s1/break_threshold)^3;
    % v2 = (s2/break_threshold)^3;

    num_points = 200;
v1 = linspace(0.1 , 10, num_points);    % 避免 T_i <= Delta_t
v2 = linspace(0.1 , 10, num_points);% 避免 tau_p <= Delta_t
[v1_grid, v2_grid] = meshgrid(v1, v2);

% 自定义聚合概率函数 


    %%%%%%% 实验参数 %%%%%%%
    Re_tau = 180;
    delta = 0.5e-2;
    nu_c = 1.8e-6;
    rho_c = 1000;
    rho_d = 866;
    sigma = 0.004;
    %%%%%%% 实验参数 %%%%%%%

    %%%%%%% 模型参数 %%%%%%%
        K_1 = 3e-4;
        K_2 = 3e-4;
        K_3 = 5e-3;
    %%%%%%% 模型参数 %%%%%%%
    
    % for i = 1:num_points
    %     v = v1(i);
    %     dissipation = Re_tau^3*nu_c^3/delta^4;
    %     breakup_rate(i) = v^(-2/9)*dissipation^(1/3)*exp(-K_3*sigma/(rho_d*dissipation^(2/3)*v^(5/9)));
    % end
    % plot(v1,breakup_rate);

    % 基础概率与尺寸乘积成正比，保证小液滴聚合概率低
    % p = min(0.0005, 0.01*(s1/breakiyou_threshold)*(s2/break_threshold));2358.mat
    
% 预分配结果矩阵
result_p = zeros(size(v1_grid));

% 逐点计算表达式值
for i = 1:num_points
    for j = 1:num_points
        dissipation = Re_tau^3*nu_c^3/delta^4;
        v = sqrt(v1_grid(i,j)*v2_grid(i,j));
        breakup_rate = v^(-2/9)*dissipation^(1/3)*exp(-K_3*sigma/(rho_d*dissipation^(2/3)*v^(5/9)));
        collision = K_1*(v1_grid(i,j).^(2/3) + v2_grid(i,j).^(2/3)).*sqrt(v1_grid(i,j).^(2/9)+v2_grid(i,j).^(2/9))*Re_tau*nu_c/(delta^(4/3));
        efficiency = exp(-K_2*nu_c^4*rho_c^2*Re_tau^3*(v1_grid(i,j).^(1/3).*v2_grid(i,j).^(1/3)./(v1_grid(i,j).^(1/3)+v2_grid(i,j).^(1/3))).^4/(sigma^2*delta^4));
        
        result_p(i,j) = collision.*efficiency/breakup_rate;
    end
end
    % figure;
    % clf;
mesh(v1_grid, v2_grid, result_p);%Coulaloglou C A, Tavlarides L L. Description of interaction processes in agitated liquid-liquid dispersions[J]. Chemical Engineering Science, 1977, 32(11): 1289-1297.
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
% set(gca, 'ZScale', 'log');
xlabel('v1')
ylabel('v2')
title(['K_1 = ',num2str(K_1),',K_2 = ',num2str(K_2)]);
colorbar
hold on



% clear,clc;
% A = rand(20000,20000);
% n = 50000;
% top_values = maxk(A(:), n);
% nth_max = top_values(n);

%%
f_v = linspace(0,1,1000);
alpha = 2;


beta  = 20;
beta1 = (betapdf(f_v, alpha,beta)+ betapdf(1-f_v, alpha,beta))/2+normpdf(f_v, 0.5, 1);
plot(f_v,beta1,DisplayName=['beta=',num2str(beta)])
hold on
% 
% beta = 10;
% beta2 = (betapdf(f_v, alpha,beta)+ betapdf(1-f_v, alpha,beta))/2;
% plot(f_v,beta2,DisplayName=['beta=',num2str(beta)])
% 
% legend();
% 
% %%
% 
% clear,clc,close all;
% % 导入数据
% re_tau = 180;
% filename = ['../data_',num2str(re_tau),'/Re_tau = ',num2str(re_tau),'.mat'];
% load(filename)
% 
% re_5210.di_less_dia = [0.19847	0.39978	0.59622	0.80178	1.00309	1.2044	1.4051	1.60641	1.80286	2.00842	2.2	2.40131	2.60201	2.80332];
% re_5210.pdf = [0.03176	0.43487	1.08767	1.00889	0.72896	0.51658	0.48737	0.24535	0.18075	0.12084	0.05521	0.03559	0.0234	0.02726];
% 
% re_7810.di_less_dia = [0.50013	0.60108	0.7063	0.80178	0.89788	1.00309	1.09858	1.2044	1.29989	1.4051	1.5012	1.60641	1.7019	1.80711	1.90321	2.00356	2.10451	2.2];
% re_7810.pdf = [0.9895	1.21899	1.02616	0.8493	1.00889	0.8493	0.71496	0.80128	0.60186	0.54754	0.33304	0.30817	0.18744	0.13809	0.09809	0.11209	0.06952	0.04147];
% 
% % 参数设置
% initial_size = 20*data.critical_value(1)*ones(1,1);     % 初始液滴大小
% break_threshold = data.critical_value(1);  % 破碎尺寸阈值
% num_steps = 2000;        % 仿真步数
% output_step = round(num_steps/100);       % 结果输出间隔
% 
% 
% load("2143.mat")
% 
% 
% final_sizes = history{end};
% figure;
% subplot(2,1,1);
% normalized_final_sizes = final_sizes/mean(final_sizes);
% h = histogram(normalized_final_sizes,"Normalization","pdf");%, 'BinWidth', 5);
% hold on
% line([break_threshold/mean(final_sizes) break_threshold/mean(final_sizes)],[0.01 1],'linewidth',2);
% plot(re_5210.di_less_dia,re_5210.pdf,'rx',"LineWidth",2,"DisplayName",'Re = 5210 (Exp.)')
% hold on
% plot(re_7810.di_less_dia,re_7810.pdf,'bx',"LineWidth",2,"DisplayName",'Re = 7810 (Exp.)')
% hold off
% set(gca, 'YScale', 'log');
% title(['最终液滴尺寸分布，均值为',num2str(mean(final_sizes))]);
% xlabel('液滴尺寸');
% ylabel('数量');
% 
% 
% subplot(2,1,2);
% num_droplets = cellfun(@length, history);
% for i = 1:length(num_droplets)
%     history_mean(i) = mean(num_droplets(1:i));
% end
% plot(1:length(num_droplets), num_droplets, '-o');
% hold on
% plot(1:length(num_droplets), history_mean, '-r',LineWidth=1);
% title('液滴数量演化');
% xlabel('时间步');
% ylabel('液滴总数');
% grid on;
