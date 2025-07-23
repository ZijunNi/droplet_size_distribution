% 本代码用于对stat文件夹中保存的各项统计数据进行合并绘图
clear,clc
stat_dir = 'stat';

re_tau = 180;
ratio = [1,2,5,10];%百分位数

re_5210.di_less_dia = [0.19847	0.39978	0.59622	0.80178	1.00309	1.2044	1.4051	1.60641	1.80286	2.00842	2.2	2.40131	2.60201	2.80332];
re_5210.pdf = [0.03176	0.43487	1.08767	1.00889	0.72896	0.51658	0.48737	0.24535	0.18075	0.12084	0.05521	0.03559	0.0234	0.02726];

re_7810.di_less_dia = [0.50013	0.60108	0.7063	0.80178	0.89788	1.00309	1.09858	1.2044	1.29989	1.4051	1.5012	1.60641	1.7019	1.80711	1.90321	2.00356	2.10451	2.2];
re_7810.pdf = [0.9895	1.21899	1.02616	0.8493	1.00889	0.8493	0.71496	0.80128	0.60186	0.54754	0.33304	0.30817	0.18744	0.13809	0.09809	0.11209	0.06952	0.04147];


figure;
    plot(re_5210.di_less_dia,re_5210.pdf,'rs',"LineWidth",2,"DisplayName",'Re = 5210 (Exp.)')
    hold on
    plot(re_7810.di_less_dia,re_7810.pdf,'bx',"LineWidth",2,"DisplayName",'Re = 7810 (Exp.)')





for i = 1:4
    filename =fullfile(stat_dir,['stat_re_tau_',num2str(re_tau),'_DNS_one_para_baseline_ratio_',num2str(ratio(i)),'.mat']);
    load(filename)
    mean_size(i) = mean(data_of_droplets);
    if(1)% 重新计算PDF
        edges = [0.02:0.05:1.04 1.05:0.1:3, inf];  % 最后一个分箱捕获所有极值
        [counts, edges] = histcounts(data_of_droplets/mean_size(i), edges,"Normalization","pdf");
        centers = edges(1:end-1) + diff(edges)/2;
    end
    plot(centers*mean_size(i)/mean_size(1), counts, 'LineWidth', 1.5,DisplayName=['Theory, eta = ',num2str(ratio(i)),'%']);
end

filename =fullfile(stat_dir,'19_Jun_2025_17_38_30.mat');
load(filename);
plot(centers,counts,'LineWidth', 1.5,DisplayName='Random Perturbation');

legend();
hold off
set(gca, 'YScale', 'log');
% set(gca, 'XScale', 'log');
legend();
% title(['Final Droplets Size Distribution with Mean Value = ',num2str(mean(data_of_droplets),'%.5e')]);
xlabel('Droplet Size $D/\langle D\rangle$','Interpreter','latex');
xlim([0 4])
ylabel('PDF');
