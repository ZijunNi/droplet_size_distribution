f_v = linspace(0,1,1000);
alpha = 2;


beta  = 100;
beta1 = (betapdf(f_v, alpha,beta)+ betapdf(1-f_v, alpha,beta))/2;
plot(f_v,beta1,DisplayName=['beta=',num2str(beta)])
hold on

beta = 10;
beta2 = (betapdf(f_v, alpha,beta)+ betapdf(1-f_v, alpha,beta))/2;
plot(f_v,beta2,DisplayName=['beta=',num2str(beta)])

legend();

%%

clear,clc,close all;
% 导入数据
re_tau = 180;
filename = ['../data_',num2str(re_tau),'/Re_tau = ',num2str(re_tau),'.mat'];
load(filename)

re_5210.di_less_dia = [0.19847	0.39978	0.59622	0.80178	1.00309	1.2044	1.4051	1.60641	1.80286	2.00842	2.2	2.40131	2.60201	2.80332];
re_5210.pdf = [0.03176	0.43487	1.08767	1.00889	0.72896	0.51658	0.48737	0.24535	0.18075	0.12084	0.05521	0.03559	0.0234	0.02726];

re_7810.di_less_dia = [0.50013	0.60108	0.7063	0.80178	0.89788	1.00309	1.09858	1.2044	1.29989	1.4051	1.5012	1.60641	1.7019	1.80711	1.90321	2.00356	2.10451	2.2];
re_7810.pdf = [0.9895	1.21899	1.02616	0.8493	1.00889	0.8493	0.71496	0.80128	0.60186	0.54754	0.33304	0.30817	0.18744	0.13809	0.09809	0.11209	0.06952	0.04147];

% 参数设置
initial_size = 20*data.critical_value(1)*ones(1,1);     % 初始液滴大小
break_threshold = data.critical_value(1);  % 破碎尺寸阈值
num_steps = 2000;        % 仿真步数
output_step = round(num_steps/100);       % 结果输出间隔


load("2143.mat")


final_sizes = history{end};
figure;
subplot(2,1,1);
normalized_final_sizes = final_sizes/mean(final_sizes);
h = histogram(normalized_final_sizes,"Normalization","pdf");%, 'BinWidth', 5);
hold on
line([break_threshold/mean(final_sizes) break_threshold/mean(final_sizes)],[0.01 1],'linewidth',2);
plot(re_5210.di_less_dia,re_5210.pdf,'rx',"LineWidth",2,"DisplayName",'Re = 5210 (Exp.)')
hold on
plot(re_7810.di_less_dia,re_7810.pdf,'bx',"LineWidth",2,"DisplayName",'Re = 7810 (Exp.)')
hold off
set(gca, 'YScale', 'log');
title(['最终液滴尺寸分布，均值为',num2str(mean(final_sizes))]);
xlabel('液滴尺寸');
ylabel('数量');


subplot(2,1,2);
num_droplets = cellfun(@length, history);
for i = 1:length(num_droplets)
    history_mean(i) = mean(num_droplets(1:i));
end
plot(1:length(num_droplets), num_droplets, '-o');
hold on
plot(1:length(num_droplets), history_mean, '-r',LineWidth=1);
title('液滴数量演化');
xlabel('时间步');
ylabel('液滴总数');
grid on;
