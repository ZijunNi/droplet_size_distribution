clear,clc,close all;
% 导入数据
re_tau = 180;
filename = ['../data_',num2str(re_tau),'/Re_tau = ',num2str(re_tau),'.mat'];
load(filename)

re_5210.di_less_dia = [0.19847	0.39978	0.59622	0.80178	1.00309	1.2044	1.4051	1.60641	1.80286	2.00842	2.2	2.40131	2.60201	2.80332];
re_5210.pdf = [0.03176	0.43487	1.08767	1.00889	0.72896	0.51658	0.48737	0.24535	0.18075	0.12084	0.05521	0.03559	0.0234	0.02726];

re_7810.di_less_dia = [0.50013	0.60108	0.7063	0.80178	0.89788	1.00309	1.09858	1.2044	1.29989	1.4051	1.5012	1.60641	1.7019	1.80711	1.90321	2.00356	2.10451	2.2];
re_7810.pdf = [0.9895	1.21899	1.02616	0.8493	1.00889	0.8493	0.71496	0.80128	0.60186	0.54754	0.33304	0.30817	0.18744	0.13809	0.09809	0.11209	0.06952	0.04147];


%%%%%%%%%% 可调参数：初始状态和程序设置 %%%%%%%%%%
initial_size = 10*data.critical_value(1)*rand(20,1);% 初始液滴大小分布
num_steps = 2000;                                   % 步数
output_step = round(num_steps/100);                 % 结果输出间隔
%%%%%%%%%% 可调参数：初始状态和程序设置 %%%%%%%%%%

% 确保离散相总体积不变
fixed_total_volume = (20*data.critical_value(1))^(3);% 固定总体积为8000*V_crit
total_volume = sum(initial_size.^3);
normalized_size = initial_size./(total_volume^(1/3));
initial_size = normalized_size*(fixed_total_volume)^(1/3);
% 已验证：最终平均直径与离散相总体积有关




break_threshold = data.critical_value(1);  % 破碎尺寸阈值
% 初始化液滴系统
droplets = initial_size;
history = cell(num_steps+1, 1);
history{1} = droplets;

% 主循环
num_met = 0;
for step = 1:num_steps
    % 破碎过程
    new_droplets = [];
    for i = 1:length(droplets)
        s = droplets(i);
        if s > break_threshold

            %%%%%%%%%% 可调参数：分裂比例模型 %%%%%%%%%%
            a = 2;
            b = 10*s/break_threshold;%ratio(2);    
            %%%%%%%%%% 可调参数：分裂比例模型 %%%%%%%%%%

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
    
    % 考虑所有可能配对
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
                %%%%%%%%%% 可调参数：聚合概率模型 %%%%%%%%%%
                prob = aggregation_prob(s1, s2, break_threshold);
                %%%%%%%%%% 可调参数：聚合概率模型 %%%%%%%%%%
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
    num_droplets(step) = length(droplets);
    history_mean(step) = mean(num_droplets);

    history_mean_size(step) = mean(droplets);
    history_mean_std(step) = std(droplets);

    % 退出条件：连续多步以上标准差小于均值的百分之一
    if(step>output_step)
        error = std(num_droplets(end-output_step:end));
        criteria = 0.01*mean(num_droplets(end-output_step:end));
        
        if (error < criteria&&num_met>floor(num_steps/100))
            break; 
        elseif(error < criteria)
            num_met = num_met + 1;
        elseif(error > criteria)
            num_met = 0;
        end
    else
        error = Inf;
    end


    % 阶段性输出
    if mod(step, output_step) == 0
        fprintf('Step %d: %d droplets, Error: %.5e\n',...
                step, length(droplets),error);
    end
end

%% 绘图
final_sizes = droplets;
figure;
subplot(3,1,1);
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


subplot(3,1,2);

plot(1:length(num_droplets), num_droplets, '-o');
% hold on
% plot(1:length(num_droplets), history_mean, '-r',LineWidth=1);
% title('液滴数量演化');
xlabel('Steps');
ylabel('Numbers of Droplets');
grid on;

subplot(3,1,3);
semilogy(1:length(num_droplets), history_mean_size, '-x',DisplayName='Mean Droplets Size');
hold on
semilogy(1:length(num_droplets), history_mean_std, '-o',DisplayName='Std. Deviation of Droplets Size');
hold off
legend();
title('Statisitic Values');
xlabel('Steps');
ylabel('Statistic Values');
grid on;



% 自定义聚合概率函数 
function p = aggregation_prob(s1, s2, break_threshold)
    v1 = (s1/break_threshold)^3;
    v2 = (s2/break_threshold)^3;
    % 基础概率与尺寸乘积成正比，保证小液滴聚合概率低
    % p = min(0.0005, 0.01*(s1/breakiyou_threshold)*(s2/break_threshold));  % 2358.mat
    p = 0.0001*(v1^(2/3) + v2^(2/3))*sqrt(v1^(2/9)+v2^(2/9));%Coulaloglou C A, Tavlarides L L. Description of interaction processes in agitated liquid-liquid dispersions[J]. Chemical Engineering Science, 1977, 32(11): 1289-1297.

end
