% 本代码处理多参数实验数据，若要处理多次重复实验的数据，请使用coalescence_breakup_data_process.m

clc,clear

% 导入数据
re_tau = 180;
filename = ['../data_',num2str(re_tau),'/Re_tau = ',num2str(re_tau),'.mat'];
load(filename)

% res{1} = load("./data_re_tau_180/data_re_tau=180_DNS_full_model.mat","history").history;
% res{2} = load("./data_re_tau_180/data_re_tau=180_DNS_no_break_rate.mat","history").history;
% res{3} = load("./data_re_tau_180/data_re_tau=180_DNS_no_break_rate_and_efficiency.mat","history").history;
% res{4} = load("./data_re_tau_180/data_re_tau=180_DNS_no_break_rate_and_efficiency_and_coal_rate.mat","history").history;
res{1} = load("./data_re_tau_180/data_re_tau=180_DNS_vol_frac_1.mat","history").history;
res{2} = load("./data_re_tau_180/data_re_tau=180_DNS_vol_frac_5.mat","history").history;
res{3} = load("./data_re_tau_180/data_re_tau=180_DNS_vol_frac_10.mat","history").history;
res{4} = load("./data_re_tau_180/data_re_tau=180_DNS_4_vol_frac_2.mat","history").history;


ratio = data.ratio;

re_5210.di_less_dia = [0.19847	0.39978	0.59622	0.80178	1.00309	1.2044	1.4051	1.60641	1.80286	2.00842	2.2	2.40131	2.60201	2.80332];
re_5210.pdf = [0.03176	0.43487	1.08767	1.00889	0.72896	0.51658	0.48737	0.24535	0.18075	0.12084	0.05521	0.03559	0.0234	0.02726];

re_7810.di_less_dia = [0.50013	0.60108	0.7063	0.80178	0.89788	1.00309	1.09858	1.2044	1.29989	1.4051	1.5012	1.60641	1.7019	1.80711	1.90321	2.00356	2.10451	2.2];
re_7810.pdf = [0.9895	1.21899	1.02616	0.8493	1.00889	0.8493	0.71496	0.80128	0.60186	0.54754	0.33304	0.30817	0.18744	0.13809	0.09809	0.11209	0.06952	0.04147];

%%

figure;
    plot(re_5210.di_less_dia,re_5210.pdf,'rs',"LineWidth",2,"DisplayName",'Re = 5210 (Exp.)')
    hold on
    plot(re_7810.di_less_dia,re_7810.pdf,'bx',"LineWidth",2,"DisplayName",'Re = 7810 (Exp.)')
    
for index_ratio = 1:4

    final_sizes = cell2mat(res{index_ratio}(end-100:10:end)');%合并最后几步的液滴群获得平滑pdf
    normalized_final_sizes = final_sizes/mean(final_sizes);

    if(1)%处理离群值较多之情况
        % data = normalized_final_sizes;

        edges = [0.02:0.1:1 1.05:0.1:3, inf];  % 最后一个分箱捕获所有极值
        % [counts, edges] = histcounts(data, edges);
 
        [counts, edges] = histcounts(normalized_final_sizes, edges,"Normalization","pdf");
    else
        [counts, edges] = histcounts(normalized_final_sizes,75,"Normalization","pdf");
    end

    centers = edges(1:end-1) + diff(edges)/2;
    plot(centers, counts, 'LineWidth', 1.5,DisplayName=['Theory, \eta = ',num2str(100*ratio(index_ratio)),'%']);
    hold on
    % line([break_threshold/mean(final_sizes) break_threshold/mean(final_sizes)],[0.01 1],'linewidth',2,Displayname='Critical Value');
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    legend();
    title(['Final Droplets Size Distribution with Mean Value = ',num2str(mean(final_sizes),'%.5e')]);
    xlabel('Droplet Size $D/\langle D\rangle$','Interpreter','latex');
    xlim([0 4])
    ylabel('PDF');

end

plot([0.5,1.8,NaN],[0.5,1.8,NaN].^(-3/2)*0.65,DisplayName='k = -3/2',LineWidth=3)
plot([1.5,3],[1.5,3].^(-10/3)*0.6,DisplayName='k = -10/3')
hold off

%%

vol_frac = [1,5,10,20];
for index_ratio = 1:3
mean_val(index_ratio) = mean(cell2mat(res{index_ratio}(end-100:5:end)'));

end
mean_val(4) = 15.649*4e-3/22.208+0.24e-1;

figure
plot(vol_frac,mean_val*1e+4,'x',MarkerSize=10,LineWidth=2);
ylim([0 400])
xlim([0 25])
ylabel('$\langle D \rangle/\rm{\mu m}$','Interpreter','latex')
xlabel('$\rm{Volume\; Fraction}\;\phi/\%$','Interpreter','latex')
% hold on 
% semilogx(vol_frac,log(vol_frac)*0.0039+mean_val(1))
% ylim([0.014 0.015]);
% xlim([0 0.12]);


%%
% 参数设置
initial_size = 20*data.critical_value(1)*ones(1,1);     % 初始液滴大小
break_threshold = data.critical_value(1);  % 破碎尺寸阈值
num_steps = 2000;        % 步数
output_step = round(num_steps/100);       % 结果输出间隔






% 定义示例元胞数组（可根据实际情况替换）
cellArray = history(end);


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