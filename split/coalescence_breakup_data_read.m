clc,clear

% 导入数据
re_tau = 180;
filename = ['../data_',num2str(re_tau),'/Re_tau = ',num2str(re_tau),'.mat'];
load(filename)
load("coal_split_1.mat");

re_5210.di_less_dia = [0.19847	0.39978	0.59622	0.80178	1.00309	1.2044	1.4051	1.60641	1.80286	2.00842	2.2	2.40131	2.60201	2.80332];
re_5210.pdf = [0.03176	0.43487	1.08767	1.00889	0.72896	0.51658	0.48737	0.24535	0.18075	0.12084	0.05521	0.03559	0.0234	0.02726];

re_7810.di_less_dia = [0.50013	0.60108	0.7063	0.80178	0.89788	1.00309	1.09858	1.2044	1.29989	1.4051	1.5012	1.60641	1.7019	1.80711	1.90321	2.00356	2.10451	2.2];
re_7810.pdf = [0.9895	1.21899	1.02616	0.8493	1.00889	0.8493	0.71496	0.80128	0.60186	0.54754	0.33304	0.30817	0.18744	0.13809	0.09809	0.11209	0.06952	0.04147];


% 参数设置
initial_size = 20*data.critical_value(1)*ones(1,1);     % 初始液滴大小
break_threshold = data.critical_value(1);  % 破碎尺寸阈值
num_steps = 2000;        % 步数
output_step = round(num_steps/100);       % 结果输出间隔




history(cellfun(@isempty,history))=[];

% 定义示例元胞数组（可根据实际情况替换）
cellArray = history(end-20:end);


% 将每个元胞元素转换为列向量并合并
mergedVector = cell2mat(cellfun(@(x) x(:), cellArray, 'UniformOutput', false));

% 显示结果
final_sizes = mergedVector;
figure;

normalized_final_sizes = final_sizes/mean(final_sizes);
h = histogram(normalized_final_sizes,"Normalization","pdf",DisplayName='Theory');
hold on
line([break_threshold/mean(final_sizes) break_threshold/mean(final_sizes)],[0.01 1],'linewidth',2,Displayname='Critical Value');
plot(re_5210.di_less_dia,re_5210.pdf,'rx',"LineWidth",2,"DisplayName",'Re = 5210 (Exp.)')
hold on
plot(re_7810.di_less_dia,re_7810.pdf,'bx',"LineWidth",2,"DisplayName",'Re = 7810 (Exp.)')
hold off
set(gca, 'YScale', 'log');
legend();
title(['Final Droplets Size Distribution with Mean Value = ',num2str(mean(final_sizes))]);
xlabel('Droplet Size $D/\langle D\rangle$','Interpreter','latex');
ylabel('PDF');