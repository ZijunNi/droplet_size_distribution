function [edges,counts] = whole_model(b,K_1)


%%%%%%%%%% 可调参数：初始状态和程序设置 %%%%%%%%%%
initial_size = 5 + rand(200,1);             % 初始液滴大小分布
num_steps = 200;                            % 步数
output_step = round(num_steps/100);         % 结果输出间隔
%%%%%%%%%% 可调参数：初始状态和程序设置 %%%%%%%%%%

    vol_frac = [1,5,10,20];
    % 导入数据
    re_tau = 180;
    % filename = ['../data_',num2str(re_tau),'/Re_tau = ',num2str(re_tau),'.mat'];
    % load(filename)

    data.critical_value(1) = 0.0141276673007496;

    clear history num_droplets history_mean history_mean_size history_mean_std;

    index_ratio = 1;%固定百分位数
    index_vol_frac = 1;% 固定体积分数

    fixed_total_volume = vol_frac(index_vol_frac)*(15*data.critical_value(1))^3;%0.01*pi()*(3.5^2-2.5^2)*7.5;% 固定总体积为0.01*总体积，总体积来自Sun Chao实验中的装置总体积
    total_volume = sum(initial_size.^3);
    normalized_size = initial_size./(total_volume^(1/3));
    initial_size = normalized_size*(fixed_total_volume)^(1/3);% 已验证：最终平均直径与离散相总体积有关，但与初始分布无关

    close all;
    break_threshold = data.critical_value(index_ratio);  % 破碎尺寸阈值
    
    % 初始化液滴系统
    droplets = initial_size;
    history = cell(num_steps+1, 1);
    history{1} = droplets;
    
    % 主循环
    num_met = 100;%达到目标后停止的步数
    
    
    flag_met = false;
    target_step = Inf;
    extra_steps = 100;  % 设置需要的额外运行步数
    
    step = 1;
    while(step <= num_steps)
        % 破碎过程
        new_droplets = [];
        for i = 1:length(droplets)
            s = droplets(i);
            % break_threshold = data.critical_value(index_ratio)*(1-rand(1,1)*0.01);% 引入临界值的随机漂移
            if s > break_threshold
    
                %%%%%%%%%% 可调参数：分裂比例模型 %%%%%%%%%%
                a = 2;
                % b = 9;%*s/break_threshold;%20;%10*s/break_threshold;%ratio(2);  %20常数效果也不错
                
                %%%%%%%%%% 可调参数：分裂比例模型 %%%%%%%%%%
    
                ratio_val = custom_beta_rnd(a,b,0.1);%custom_rand_ao_shape(1,0.15,0.2);%custom_beta_rnd(a,b,0.1);%敏感性验证：需要较悬殊比例分离才可能较吻合结果
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
                    [prob, K] = model_aggregation_prob(s1, s2, break_threshold,re_tau,K_1);
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
        if(step > output_step&&num_droplets(step)>40)
            error = std(num_droplets(end-output_step:end));
            criteria = 0.01*mean(num_droplets(end-output_step:end));
            
            if (error < criteria && num_met > floor(num_steps/100))
                if ~flag_met  % 首次满足条件时设置目标步数
                    flag_met = true;
                    target_step = step + extra_steps;
                end
            elseif(error < criteria)
                num_met = num_met + 1;
            elseif(error > criteria)
                num_met = 0;
            end

        else
            error = Inf;
        end
        
        % 在每次循环迭代最后添加退出检查
        if flag_met && (step >= target_step)
            break; 
        end
    
        if(length(droplets)<5)
            warning('液滴太少，重新迭代。');
            pause(0.1);
            step = 1;
            continue;
        end
    
    
        % 阶段性输出
        if mod(step, output_step) == 0
            fprintf('Step %d: %d droplets, Error: %.5e, Mean Size: %.5e\n',...
                    step, length(droplets),error,mean(droplets));
        end
        step = step + 1;
    end
    
    history(cellfun(@isempty,history))=[];
    
    % 绘图
        if(1)%绘图开关
            % figure;
            [edges,counts] = plot_droplets_stat(history,20,5,1);
        end
    % 数据保存
        date_str = datestr(now);
        % 需要替换的字符列表
        chars_to_replace = [':', ' ','-'];
        replace_positions = ismember(date_str, chars_to_replace);
        % 执行替换
        date_str(replace_positions) = '_';
        save(['./data_re_tau_',num2str(re_tau),'/data_re_tau_',num2str(re_tau),'_',date_str,'.mat'],"history","K","b");
    

    function [p,K] = model_aggregation_prob(s1, s2, break_threshold, Re_tau,K_1)
    
        %%%%%%% 实验参数 %%%%%%%
            delta = 0.5e-2;
            nu_c = 1.8e-6;
            rho_c = 1000;
            rho_d = 866;
            sigma = 0.004;
        %%%%%%% 实验参数 %%%%%%%
    
        %%%%%%% 模型参数 %%%%%%%
            % K_1 = 3.4e-4;
            K_2 = 0;%3e-4;
            K_3 = 0;%5e-3;
            K = [K_1,K_2,K_3];
        %%%%%%% 模型参数 %%%%%%%
    
        v1 = (s1/break_threshold)^3;
        v2 = (s2/break_threshold)^3;
        % 基础概率与尺寸乘积成正比，保证小液滴聚合概率低
        % p = min(0.0005, 0.01*(s1/breakiyou_threshold)*(s2/break_threshold));2358.mat
        dissipation = Re_tau^3*nu_c^3/delta^4;
    
    
        collision = K_1*(v1^(2/3) + v2^(2/3))*sqrt(v1^(2/9)+v2^(2/9))*dissipation^(1/3);%取为常数1e-4（K2K3均为0）时效果较好
        efficiency = exp(-K_2*nu_c*rho_c^2*dissipation*(v1^(1/3)*v2^(1/3)/(v1^(1/3)+v2^(1/3)))^4/(sigma^2));
        
        v = sqrt(v1*v2);
        breakup_rate = v^(-2/9)*dissipation^(1/3)*exp(-K_3*sigma/(rho_d*dissipation^(2/3)*v^(5/9)));
    
        p = collision*efficiency/breakup_rate;
        if(p>0.99)
            warning('最大概率接近1，无法形成液滴群。')
        end
        %Coulaloglou C A, Tavlarides L L. Description of interaction processes in agitated liquid-liquid dispersions[J]. Chemical Engineering Science, 1977, 32(11): 1289-1297.
    
    end


end