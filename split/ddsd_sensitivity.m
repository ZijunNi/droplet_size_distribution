clear,clc



M = 20;S = 1;
edges = [0.02:0.1:1.04, 1.05:0.1:3, inf];

% Experiment
re_5210.di_less_dia = [0.19847	0.39978	0.59622	0.80178	1.00309	1.2044	1.4051	1.60641	1.80286	2.00842	2.2	2.40131	2.60201	2.80332];
re_5210.pdf = [0.03176	0.43487	1.08767	1.00889	0.72896	0.51658	0.48737	0.24535	0.18075	0.12084	0.05521	0.03559	0.0234	0.02726];

re_7810.di_less_dia = [0.50013	0.60108	0.7063	0.80178	0.89788	1.00309	1.09858	1.2044	1.29989	1.4051	1.5012	1.60641	1.7019	1.80711	1.90321	2.00356	2.10451	2.2];
re_7810.pdf = [0.9895	1.21899	1.02616	0.8493	1.00889	0.8493	0.71496	0.80128	0.60186	0.54754	0.33304	0.30817	0.18744	0.13809	0.09809	0.11209	0.06952	0.04147];

plot(re_5210.di_less_dia, re_5210.pdf, 'ro', "LineWidth", 2, "DisplayName", 'Re = 5210 (Exp.)');hold on;
plot(re_7810.di_less_dia, re_7810.pdf, 'bx', "LineWidth", 2, "DisplayName", 'Re = 7810 (Exp.)');



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
edges = [0.02:0.1:0.84,0.85,1,1.19, 1.2:0.1:3, inf];
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
 