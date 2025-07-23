clear,clc,close all;
% 导入数据
re_tau = 1000;
filename = ['../data_',num2str(re_tau),'/Re_tau = ',num2str(re_tau),'.mat'];
load(filename)

re_5210.di_less_dia = [0.19847	0.39978	0.59622	0.80178	1.00309	1.2044	1.4051	1.60641	1.80286	2.00842	2.2	2.40131	2.60201	2.80332];
re_5210.pdf = [0.03176	0.43487	1.08767	1.00889	0.72896	0.51658	0.48737	0.24535	0.18075	0.12084	0.05521	0.03559	0.0234	0.02726];

re_7810.di_less_dia = [0.50013	0.60108	0.7063	0.80178	0.89788	1.00309	1.09858	1.2044	1.29989	1.4051	1.5012	1.60641	1.7019	1.80711	1.90321	2.00356	2.10451	2.2];
re_7810.pdf = [0.9895	1.21899	1.02616	0.8493	1.00889	0.8493	0.71496	0.80128	0.60186	0.54754	0.33304	0.30817	0.18744	0.13809	0.09809	0.11209	0.06952	0.04147];


%%%%%%%%%% 可调参数：初始状态和程序设置 %%%%%%%%%%
initial_size = 1+rand(1,400);                                                    % 初始液滴大小分布
num_steps = 500;                                               % 最大步数
output_step = round(num_steps/100);            % 结果输出间隔
%%%%%%%%%% 可调参数：初始状态和程序设置 %%%%%%%%%%


index_eta = 1;

% for index_eta = 1:4

% 确保离散相总体积不变
% phi = 0.01;
% total_volume_of_domain = pi()*(3.5^2-2.5^2)*7.5;% Sun Chao 实验中的装置总体积
phi = 1;
fixed_total_volume = phi*0.03^3;%0.01*pi()*(3.5^2-2.5^2)*7.5;% 固定总体积为0.01*总体积，总体积来自Sun Chao实验中的装置总体积
total_volume = sum(initial_size.^3);
normalized_size = initial_size./(total_volume^(1/3));
initial_size = normalized_size*(fixed_total_volume)^(1/3);% 已验证：最终平均直径与离散相总体积有关

initial_size = initial_size/data.critical_value(index_eta);% 初始尺寸归一化

clear history_mean_std history history_mean_size history_mean droplets

break_threshold = 1;%data.critical_value(index_eta);  % 破碎尺寸阈值，归一化后为1
% 初始化液滴系统
droplets = initial_size;
history = cell(num_steps+1, 1);
history{1} = droplets;

tic
% 主循环
num_met = 0;
steps = 1;
break_record = [];
while (steps < num_steps)
    % 破碎过程
    new_droplets = [];
    for i = 1:length(droplets)
        s = droplets(i);

        if s > break_threshold
            % break_record = [break_record, s];
            %%%%%%%%%% 可调参数：分裂比例模型 %%%%%%%%%%
            a = 2;
            % b = 10*s;
            % ratio_val = custom_beta_rnd(a,b,1);
            ratio_val = calculate_ratio_value(s/break_threshold,re_tau);%custom_beta_rnd(a,b,1);
            %%%%%%%%%% 可调参数：分裂比例模型 %%%%%%%%%%
            
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
                [prob, K] = aggregation_prob(s1, s2, data.critical_value(index_eta), re_tau);
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

    if(length(droplets)<100)% 液滴群未破碎完全
        warning('液滴群规模过小，重新执行破碎。')
        steps = 1;
        continue;
    end
    
    % 记录当前状态
    history{steps+1} = droplets;
    num_droplets(steps) = length(droplets);
    history_mean(steps) = mean(num_droplets);

    history_mean_size(steps) = mean(droplets);
    history_mean_std(steps) = std(droplets);

    % 退出条件：连续100步以上标准差小于均值的百分之一
    if(steps>output_step)
        error = std(num_droplets(end-output_step:end));
        criteria = 0.01*mean(num_droplets(end-output_step:end));
        
        if (error < criteria&&num_met>50)% 退出条件：连续100步以上标准差小于均值的百分之一
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
    if mod(steps, output_step) == 0
        fprintf('Step %d: %d droplets, Error: %.5e, Mean value: %.5e\n',...
                steps, length(droplets),error,mean(droplets));
    end
    steps = steps+1;
end
toc
history(cellfun(@isempty,history))=[];

%% 绘图

if(1) 
    M = 30;
    S = 1; 
    plot_droplets_stat(history,M,S,1);
    drawnow;
    saveas(gca, 'result.jpg');
    
    if(1)
        mailme('聚并参数调整',{num2str(K),['Re_tau = ',num2str(re_tau)],['Volnume Fraction =',num2str(phi),'%']},'./result.jpg')
    end

end

%% 保存

% 定义需要保存的变量名称（根据实际情况修改）
targetVars = {'a', 'history', 'K','index_eta','re_tau','data','phi'}; % 替换为你的变量名

% 生成格式化的时间字符串（安全文件名格式：YYYY-MM-DD_HH-MM-SS）
timeStr = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM-SS');
fileName = ['auto_saved_data_', timeStr, '.mat'];

savePath = ['./data_re_tau_',num2str(re_tau)];  % 修改为你的目标路径
filePath = fullfile(savePath, fileName);
save(filePath, targetVars{:});
% 显示保存信息
fprintf('变量已保存至文件: %s\n', fileName);
disp('保存变量列表:');
disp(targetVars);

% end




