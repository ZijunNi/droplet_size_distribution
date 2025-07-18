% 本代码处理多次重复实验的数据，若要处理多参数实验数据，请使用coalescence_breakup_data_read.m

clc,clear

% 导入独立实验的结果并合并
data_of_droplets = [];
re_tau = 180;
ratio = 10;%百分位数
for i = 1:10
    %data_re_tau_180_DNS_one_para_baseline_10_ratio_1
    filename = ['./data_re_tau_',num2str(re_tau),'/data_re_tau_',num2str(re_tau),'_DNS_one_para_baseline_',num2str(i),'_ratio_',num2str(ratio),'.mat'];
    a = load(filename).history;
    % cell2mat(res{index_ratio}(end)');
    data_of_droplets = [data_of_droplets a{end-100:10:end}];
end


re_5210.di_less_dia = [0.19847	0.39978	0.59622	0.80178	1.00309	1.2044	1.4051	1.60641	1.80286	2.00842	2.2	2.40131	2.60201	2.80332];
re_5210.pdf = [0.03176	0.43487	1.08767	1.00889	0.72896	0.51658	0.48737	0.24535	0.18075	0.12084	0.05521	0.03559	0.0234	0.02726];

re_7810.di_less_dia = [0.50013	0.60108	0.7063	0.80178	0.89788	1.00309	1.09858	1.2044	1.29989	1.4051	1.5012	1.60641	1.7019	1.80711	1.90321	2.00356	2.10451	2.2];
re_7810.pdf = [0.9895	1.21899	1.02616	0.8493	1.00889	0.8493	0.71496	0.80128	0.60186	0.54754	0.33304	0.30817	0.18744	0.13809	0.09809	0.11209	0.06952	0.04147];


figure;
    plot(re_5210.di_less_dia,re_5210.pdf,'rs',"LineWidth",2,"DisplayName",'Re = 5210 (Exp.)')
    hold on
    plot(re_7810.di_less_dia,re_7810.pdf,'bx',"LineWidth",2,"DisplayName",'Re = 7810 (Exp.)')

    % load("re_tau=180_pdf.mat");
    % plot(centers, counts, 'LineWidth', 1.5,DisplayName='Theory, Re_tau = 180');
    % load("re_tau=200_pdf.mat");
    % plot(centers, counts, 'LineWidth', 1.5,DisplayName='Theory, Re_tau = 200');
    % set(gca, 'YScale', 'log');
    % % set(gca, 'XScale', 'log');
    % legend();
    % title(['Final Droplets Size Distribution with Mean Value = ',num2str(mean(final_sizes),'%.5e')]);
    % xlabel('Droplet Size $D/\langle D\rangle$','Interpreter','latex');
    % xlim([0 4])
    % ylabel('PDF');


    
    final_sizes = data_of_droplets;%合并最后几步的液滴群获得平滑pdf
    normalized_final_sizes = final_sizes/mean(final_sizes);

    if(1)%处理离群值较多之情况
        % data = normalized_final_sizes;

        edges = [0.02:0.1:1.04 1.05:0.1:3, inf];  % 最后一个分箱捕获所有极值
        % [counts, edges] = histcounts(data, edges);
 
        [counts, edges] = histcounts(normalized_final_sizes, edges,"Normalization","pdf");
    else
        [counts, edges] = histcounts(normalized_final_sizes,75,"Normalization","pdf");
    end

    % [counts, edges] = histcounts(normalized_final_sizes,"Normalization","pdf");
    centers = edges(1:end-1) + diff(edges)/2;
    plot(centers, counts, 'LineWidth', 1.5,DisplayName=['Theory']);
    ref_points = linspace(0.5,3,100);
    % plot(ref_points,(ref_points).^(-3/2),'-.',DisplayName='Slope = -3/2');
    % ref_points = linspace(2.00842,3,100);
    % plot(ref_points,(ref_points).^(-10/3),'-',DisplayName='Slope = -10/3');
    hold off
    % line([break_threshold/mean(final_sizes) break_threshold/mean(final_sizes)],[0.01 1],'linewidth',2,Displayname='Critical Value');
    set(gca, 'YScale', 'log');
    % set(gca, 'XScale', 'log');
    legend();
    title(['Final Droplets Size Distribution with Mean Value = ',num2str(mean(final_sizes),'%.5e')]);
    xlabel('Droplet Size $D/\langle D\rangle$','Interpreter','latex');
    xlim([0 4])
    ylabel('PDF');

    % 保存数据
    save_name = fullfile('stat',['stat_re_tau_',num2str(re_tau),'_DNS_one_para_baseline_ratio_',num2str(ratio),'.mat']);
    save(save_name,"centers","counts","data_of_droplets");
