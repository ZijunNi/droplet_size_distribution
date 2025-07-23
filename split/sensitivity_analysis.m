clear,clc
%% 体积分数验证
clear,clc

M = 100;S = 1;
% edges = [0.02:0.05:0.78,0.9,1,1.19, 1.2:0.1:3, inf];
lineStyles = {':', '--', '-.'}; % 线型库，确保数量 >= n

folderPath = 'data_re_tau_180';
prefix = 'auto_saved_data_2025-07-20';
fileList = dir(fullfile(folderPath, [prefix, '*']));


% 检查是否找到文件
if isempty(fileList)
    warning('未找到以 "%s" 开头的文件', prefix);
    return;
end

% 遍历所有匹配的文件
validCount = 0;
for i = 1:length(fileList)
    % 跳过文件夹
    if fileList(i).isdir
        continue;
    end
    
    % 获取完整文件路径
    filePath = fullfile(fileList(i).folder, fileList(i).name);
    

        % 读取文件内容 (自动处理文本/二进制)
        load(filePath);
        % if fid == -1
        %     error('无法打开文件: %s', filePath);
        % end

        N = length(history);
        start_idx = max(1, N - S*(M-1));  % 防止索引小于1
        selected_history = history(start_idx:S:N);
        final_sizes = cell2mat(selected_history');
        % mean(final_sizes) = mean(final_sizes);
        normalized_final_sizes = final_sizes*1.1 / mean(final_sizes);
        mean_value(1) = mean(final_sizes);
        [counts, edges] = histcounts(normalized_final_sizes, "Normalization", "pdf");
        centers = edges(1:end-1) + diff(edges)/2;
         lineStyle = lineStyles{mod(i-1, numel(lineStyles)) + 1};
        plot(centers, counts,'LineStyle',lineStyle,'Color',[mod(0.8*i,1) 0 0] , 'LineWidth', 1.5,'DisplayName',['$\phi=',num2str(phi),'\%$']);
        hold on
        
        

end


% 设置坐标轴
set(gca, 'YScale', 'log');
legend('Location', 'northeast', 'Interpreter', 'latex');
% title(sprintf('Final Droplet Size Distribution \n Mean Size = %.5e and %6d Droplets', mean(final_sizes),length(final_sizes)));
xlabel('Normalized Droplet Size $D/\langle D\rangle$', 'Interpreter', 'latex');
ylabel('Probability Density', 'Interpreter', 'latex');
ylim([0.005 1.5]);
xlim([0 3])
% text(0,1.2,'(b)',Interpreter='latex',FontSize=20)
daspect([1/(1.5-0.005) (1/3) 1])
% 
% subplot(1,2,2)
% plot([0.01,0.05,0.1,0.2],mean_value,'o','Color',[0.8 0 0] , 'LineWidth', 1.5,'DisplayName','$\phi=20\%$');
% ylim([0.013 0.022]);
% xlim([0 0.3])
% daspect([1/(0.022-0.013) (1/0.3) 1])




%% 绘制模型完整性验证
% Experiment
re_5210.di_less_dia = [0.19847	0.39978	0.59622	0.80178	1.00309	1.2044	1.4051	1.60641	1.80286	2.00842	2.2	2.40131	2.60201	2.80332];
re_5210.pdf = [0.03176	0.43487	1.08767	1.00889	0.72896	0.51658	0.48737	0.24535	0.18075	0.12084	0.05521	0.03559	0.0234	0.02726];

re_7810.di_less_dia = [0.50013	0.60108	0.7063	0.80178	0.89788	1.00309	1.09858	1.2044	1.29989	1.4051	1.5012	1.60641	1.7019	1.80711	1.90321	2.00356	2.10451	2.2];
re_7810.pdf = [0.9895	1.21899	1.02616	0.8493	1.00889	0.8493	0.71496	0.80128	0.60186	0.54754	0.33304	0.30817	0.18744	0.13809	0.09809	0.11209	0.06952	0.04147];

plot(re_5210.di_less_dia, re_5210.pdf, 'ro', "LineWidth", 2, "DisplayName", 'Re = 5210 (Exp.)');hold on;
plot(re_7810.di_less_dia, re_7810.pdf, 'bx', "LineWidth", 2, "DisplayName", 'Re = 7810 (Exp.)');

% ratio = [1 2 5 10];


% K_2=K_3=0
M = 100;S = 2;
edges = [0.02:0.1:0.84,0.85,1,1.19, 1.2:0.1:3, inf];
file_name = 'data_re_tau_180_27_Jun_2025_11_11_21';
load(['./data_re_tau_180/',file_name,'.mat'])
N = length(history);
start_idx = max(1, N - S*(M-1));  % 防止索引小于1
selected_history = history(start_idx:S:N);
final_sizes = cell2mat(selected_history');
mean_size = mean(final_sizes);
normalized_final_sizes = final_sizes / mean(final_sizes);
[counts, edges] = histcounts(normalized_final_sizes, edges, "Normalization", "pdf");
centers = edges(1:end-1) + diff(edges)/2;
plot(centers, counts,'-.','Color',[0.3*2 0 0.8] , 'LineWidth', 1.5,'DisplayName','$Re=6642$(DNS Fields)');
hold on

M = 100;S = 1;
edges = [0.02:0.05:0.78,0.9,1,1.19, 1.2:0.1:3, inf];

file_name = 'data_re_tau_200_27_Jun_2025_22_22_44';
load(['./data_re_tau_200/',file_name,'.mat'])
N = length(history);
start_idx = max(1, N - S*(M-1));  % 防止索引小于1
selected_history = history(start_idx:S:N);
final_sizes = cell2mat(selected_history');
% mean(final_sizes) = mean(final_sizes);
normalized_final_sizes = final_sizes*1.1 / mean_size;
[counts, edges] = histcounts(normalized_final_sizes, edges, "Normalization", "pdf");
centers = edges(1:end-1) + diff(edges)/2;
plot(centers, counts,'-.','Color',[0.3 0 0] , 'LineWidth', 1.5,'DisplayName','$Re=7590$(AEH Fields)');
hold on


% 设置坐标轴
set(gca, 'YScale', 'log');
legend('Location', 'southeast', 'Interpreter', 'latex');
% title(sprintf('Final Droplet Size Distribution \n Mean Size = %.5e and %6d Droplets', mean(final_sizes),length(final_sizes)));
xlabel('Normalized Droplet Size $D/\langle D\rangle$', 'Interpreter', 'latex');
ylabel('Probability Density', 'Interpreter', 'latex');
ylim([0.005 1.5]);
xlim([0 3])
% text(0,1.2,'(b)',Interpreter='latex',FontSize=20)
daspect([1/(1.5-0.005) (1/3) 1])
 



%% 绘制模型完整性验证
edges = [0.02:0.1:0.84,0.85,1,1.19, 1.2:0.1:3, inf];
% Experiment
re_5210.di_less_dia = [0.19847	0.39978	0.59622	0.80178	1.00309	1.2044	1.4051	1.60641	1.80286	2.00842	2.2	2.40131	2.60201	2.80332];
re_5210.pdf = [0.03176	0.43487	1.08767	1.00889	0.72896	0.51658	0.48737	0.24535	0.18075	0.12084	0.05521	0.03559	0.0234	0.02726];

re_7810.di_less_dia = [0.50013	0.60108	0.7063	0.80178	0.89788	1.00309	1.09858	1.2044	1.29989	1.4051	1.5012	1.60641	1.7019	1.80711	1.90321	2.00356	2.10451	2.2];
re_7810.pdf = [0.9895	1.21899	1.02616	0.8493	1.00889	0.8493	0.71496	0.80128	0.60186	0.54754	0.33304	0.30817	0.18744	0.13809	0.09809	0.11209	0.06952	0.04147];

plot(re_5210.di_less_dia, re_5210.pdf, 'ro', "LineWidth", 2, "DisplayName", 'Re = 5210 (Exp.)');hold on;
plot(re_7810.di_less_dia, re_7810.pdf, 'bx', "LineWidth", 2, "DisplayName", 'Re = 7810 (Exp.)');

ratio = [1 2 5 10];



% K_2=K_3=0
M = 100;S = 2;
file_name = 'data_re_tau_180_27_Jun_2025_11_11_21';
load(['./data_re_tau_180/',file_name,'.mat'])
N = length(history);
start_idx = max(1, N - S*(M-1));  % 防止索引小于1
selected_history = history(start_idx:S:N);
final_sizes = cell2mat(selected_history');
normalized_final_sizes = final_sizes / mean(final_sizes);
[counts, edges] = histcounts(normalized_final_sizes, edges, "Normalization", "pdf");
centers = edges(1:end-1) + diff(edges)/2;
plot(centers, counts,'-.','Color',[0.3 0 0] , 'LineWidth', 1.5,'DisplayName','$K_2=K_3=0$');
hold on

% K_3=0
M = 20;S = 1;
file_name = 'data_re_tau_180_27_Jun_2025_12_34_53';
load(['./data_re_tau_180/',file_name,'.mat'])
N = length(history);
start_idx = max(1, N - S*(M-1));  % 防止索引小于1
selected_history = history(start_idx:S:N);
final_sizes = cell2mat(selected_history');
normalized_final_sizes = final_sizes / mean(final_sizes);
[counts, edges] = histcounts(normalized_final_sizes, edges, "Normalization", "pdf");
centers = edges(1:end-1) + diff(edges)/2;
plot(centers, counts,'-.','Color',[0.3*2 0 0.8] , 'LineWidth', 1.5,'DisplayName','$K_3=0$');
hold on

% 
M = 20;S = 1;
file_name = 'data_re_tau_180_27_Jun_2025_13_04_05';
load(['./data_re_tau_180/',file_name,'.mat'])
N = length(history);
start_idx = max(1, N - S*(M-1));  % 防止索引小于1
selected_history = history(start_idx:S:N);
final_sizes = cell2mat(selected_history');
normalized_final_sizes = final_sizes / mean(final_sizes);
[counts, edges] = histcounts(normalized_final_sizes, edges, "Normalization", "pdf");
centers = edges(1:end-1) + diff(edges)/2;
plot(centers, counts,'-.','Color',[0.3*3 0 0] , 'LineWidth', 1.5,'DisplayName','Complete Model');
hold on


% 设置坐标轴
set(gca, 'YScale', 'log');
legend('Location', 'southeast', 'Interpreter', 'latex');
% title(sprintf('Final Droplet Size Distribution \n Mean Size = %.5e and %6d Droplets', mean(final_sizes),length(final_sizes)));
xlabel('Normalized Droplet Size $D/\langle D\rangle$', 'Interpreter', 'latex');
ylabel('Probability Density', 'Interpreter', 'latex');
ylim([0.005 1.5]);
xlim([0 3])
% text(0,1.2,'(b)',Interpreter='latex',FontSize=20)
daspect([1/(1.5-0.005) (1/3) 1])
 






%% 绘制eta敏感性验证
M = 100;S = 2;
lineStyles = {':', '--', '-.','-'}; % 线型库，确保数量 >= n

edges = [0.02:0.1:0.84,0.85,1,1.19, 1.2:0.1:3, inf];
% Experiment
re_5210.di_less_dia = [0.19847	0.39978	0.59622	0.80178	1.00309	1.2044	1.4051	1.60641	1.80286	2.00842	2.2	2.40131	2.60201	2.80332];
re_5210.pdf = [0.03176	0.43487	1.08767	1.00889	0.72896	0.51658	0.48737	0.24535	0.18075	0.12084	0.05521	0.03559	0.0234	0.02726];

re_7810.di_less_dia = [0.50013	0.60108	0.7063	0.80178	0.89788	1.00309	1.09858	1.2044	1.29989	1.4051	1.5012	1.60641	1.7019	1.80711	1.90321	2.00356	2.10451	2.2];
re_7810.pdf = [0.9895	1.21899	1.02616	0.8493	1.00889	0.8493	0.71496	0.80128	0.60186	0.54754	0.33304	0.30817	0.18744	0.13809	0.09809	0.11209	0.06952	0.04147];

plot(re_5210.di_less_dia, re_5210.pdf, 'ro', "LineWidth", 2, "DisplayName", 'Re = 5210 (Exp.)');hold on;
plot(re_7810.di_less_dia, re_7810.pdf, 'bx', "LineWidth", 2, "DisplayName", 'Re = 7810 (Exp.)');

ratio = [1 2 5 10];
for i = 1:4

file_name = ['wang_ddsd_data_re_tau_180_DNS_one_para_baseline_ratio_',num2str(ratio(i))];
load(['./data_re_tau_180/',file_name,'.mat'])
N = length(history);
start_idx = max(1, N - S*(M-1));  % 防止索引小于1
selected_history = history(start_idx:S:N);
final_sizes = cell2mat(selected_history');

if(i == 1)
    mean_size = mean(final_sizes);
end
normalized_final_sizes = final_sizes / mean_size;
[counts, edges] = histcounts(normalized_final_sizes, edges, "Normalization", "pdf");
centers = edges(1:end-1) + diff(edges)/2;
 lineStyle = lineStyles{mod(i-1, numel(lineStyles)) + 1};

plot(centers, counts,'-.','Color',[0.2*i 0 0 ] ,'LineWidth', 1,'LineStyle', ':','DisplayName',['$\eta = $',num2str(ratio(i))]);

end


% 设置坐标轴
set(gca, 'YScale', 'log');
legend('Location', 'southeast', 'Interpreter', 'latex');
% title(sprintf('Final Droplet Size Distribution \n Mean Size = %.5e and %6d Droplets', mean(final_sizes),length(final_sizes)));
xlabel('Normalized Droplet Size $D/\langle D\rangle$', 'Interpreter', 'latex');
ylabel('Probability Density', 'Interpreter', 'latex');
ylim([0.005 1.5]);
xlim([0 3])
% text(0,1.2,'(b)',Interpreter='latex',FontSize=20)
daspect([1/(1.5-0.005) (1/3) 1])
 
%% 绘制DDSD敏感性验证



% Mixed Uniform Distribution
file_name = 'data_re_tau_180_27_Jun_2025_14_16_17';
load(['./data_re_tau_180/',file_name,'.mat'])
N = length(history);
start_idx = max(1, N - S*(M-1));  % 防止索引小于1
selected_history = history(start_idx:S:N);
final_sizes = cell2mat(selected_history');
normalized_final_sizes = final_sizes / mean(final_sizes);
[counts, edges] = histcounts(normalized_final_sizes, edges, "Normalization", "pdf");
centers = edges(1:end-1) + diff(edges)/2;
plot(centers, counts,'k:', 'LineWidth', 1.5,'DisplayName','Mixed Uniform Distribution');
hold on

% Bimodal DDSD
M = 20;S = 1;

file_name = 'data_re_tau_180_27_Jun_2025_17_58_57';
load(['./data_re_tau_180/',file_name,'.mat'])
N = length(history);
start_idx = max(1, N - S*(M-1));  % 防止索引小于1
selected_history = history(start_idx:S:N);
final_sizes = cell2mat(selected_history');
normalized_final_sizes = final_sizes / mean(final_sizes);
[counts, edges] = histcounts(normalized_final_sizes, edges, "Normalization", "pdf");
centers = edges(1:end-1) + diff(edges)/2;
plot(centers, counts,'k--', 'LineWidth', 1.5,'DisplayName','Bimodal DDSD');
hold on

% Bimodal DDSD
M = 100;S = 1;
file_name = 'data_re_tau_180_DNS_one_para_baseline_ratio_1';
load(['./data_re_tau_180/',file_name,'.mat'])
N = length(history);
start_idx = max(1, N - S*(M-1));  % 防止索引小于1
selected_history = history(start_idx:S:N);
final_sizes = cell2mat(selected_history');
normalized_final_sizes = final_sizes / mean(final_sizes);
[counts, edges] = histcounts(normalized_final_sizes, edges, "Normalization", "pdf");
centers = edges(1:end-1) + diff(edges)/2;
plot(centers, counts,'k-.', 'LineWidth', 1.5,'DisplayName','Bimodal DDSD with non-zero probability of equal breakup');
hold on








% 设置坐标轴
set(gca, 'YScale', 'log');
legend('Location', 'southeast', 'Interpreter', 'latex');
% title(sprintf('Final Droplet Size Distribution \n Mean Size = %.5e and %6d Droplets', mean(final_sizes),length(final_sizes)));
xlabel('Normalized Droplet Size $D/\langle D\rangle$', 'Interpreter', 'latex');
ylabel('Probability Density', 'Interpreter', 'latex');
ylim([0.005 1.5]);
xlim([0 3])
text(0,1.2,'(b)',Interpreter='latex',FontSize=20)
daspect([1/(1.5-0.005) (1/3) 1])
 