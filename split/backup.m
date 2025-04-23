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
initial_size = 10*data.critical_value(1);     % 初始液滴大小
break_threshold = data.critical_value(1);  % 破碎尺寸阈值
num_steps = 2000;        % 仿真步数
output_step = 5;       % 结果输出间隔

% 自定义聚合概率函数 
function p = aggregation_prob(s1, s2, break_threshold)
    % 基础概率与尺寸乘积成正比，保证小液滴聚合概率低
    p = min(0.8, 0.01*(s1/break_threshold)*(s2/break_threshold));  % 限制最大概率为0.8
end

% 初始化液滴系统
droplets = initial_size;
history = cell(num_steps+1, 1);
history{1} = droplets;

% 主循环
for step = 1:num_steps
    % 破碎过程
    new_droplets = [];
    for i = 1:length(droplets)
        s = droplets(i);
        if s > break_threshold
            a = 2;
            b = s;%ratio(2);    
            % y = (betapdf(x, a, b)+betapdf(1-x, a, b))/2;
            % figure;
            % plot(x, y, 'LineWidth', 2.5, ...
            %     'LineStyle', styles{mod(i-1,numel(styles))+1}, ...
            %     'Color', colors(i,:), ...
            %     'DisplayName', sprintf('α=%.1f, β=%.1f',a,b));
            
            ratio_val = custom_beta_rnd(a,b,1);
            % 采用自定义对称beta分布
            
            split_val1 = s*(1-ratio_val)^(1/3);
            split_val2 = s*ratio_val^(1/3);
            new_droplets = [new_droplets, split_val1, split_val2];
        else
            new_droplets = [new_droplets, s];
        end
    end
    droplets = new_droplets;
    
    % 增强型聚合过程 (考虑所有可能配对)
    n = length(droplets);
    if n >= 2
        % 生成所有可能配对并随机排序
        [i,j] = find(triu(ones(n),1));
        pairs = [i,j];
        pairs = pairs(randperm(size(pairs,1)),:);
        
        merged = false(1,n);
        new_droplets_agg = [];
        
        % 遍历所有可能配对
        for k = 1:size(pairs,1)
            idx1 = pairs(k,1);
            idx2 = pairs(k,2);
            
            if ~merged(idx1) && ~merged(idx2)
                s1 = droplets(idx1);
                s2 = droplets(idx2);
                
                % 计算聚合概率
                prob = aggregation_prob(s1, s2, break_threshold)
                if rand() < prob
                    % 成功聚合
                    new_droplets_agg = [new_droplets_agg, (s1^3+s2^3)^(1/3)];
                    merged(idx1) = true;
                    merged(idx2) = true;
                end
            end
        end
        
        % 添加未聚合的液滴
        droplets = [new_droplets_agg, droplets(~merged)];
    end
    
    % 记录当前状态
    history{step+1} = droplets;
    
    % 阶段性输出
    if mod(step, output_step) == 0
        fprintf('Step %d: %d droplets\n',...
                step, length(droplets));
    end
end

% 结果可视化
final_sizes = droplets;
figure;
subplot(2,1,1);
final_sizes = final_sizes/mean(final_sizes);
h = histogram(final_sizes);%, 'BinWidth', 5);
hold on
line([break_threshold/mean(final_sizes) break_threshold/mean(final_sizes)],[min(h.Values) max(h.Values)],'linewidth',2);
plot(re_5210.di_less_dia,re_5210.pdf,'rx',"LineWidth",2,"DisplayName",'Re = 5210 (Exp.)')
hold on
plot(re_7810.di_less_dia,re_7810.pdf,'bx',"LineWidth",2,"DisplayName",'Re = 7810 (Exp.)')
hold off
set(gca, 'YScale', 'log');
title('最终液滴尺寸分布');
xlabel('液滴尺寸');
ylabel('数量');

subplot(2,1,2);
num_droplets = cellfun(@length, history);
plot(0:num_steps, num_droplets, '-o');
title('液滴数量演化');
xlabel('时间步');
ylabel('液滴总数');
grid on;
