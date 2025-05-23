% 示例使用
clear,clc,clf
re_tau = 180;
filename = ['../data_',num2str(re_tau),'/Re_tau = ',num2str(re_tau),'.mat'];
load(filename)

re_5210.di_less_dia = [0.19847	0.39978	0.59622	0.80178	1.00309	1.2044	1.4051	1.60641	1.80286	2.00842	2.2	2.40131	2.60201	2.80332];
re_5210.pdf = [0.03176	0.43487	1.08767	1.00889	0.72896	0.51658	0.48737	0.24535	0.18075	0.12084	0.05521	0.03559	0.0234	0.02726];

re_7810.di_less_dia = [0.50013	0.60108	0.7063	0.80178	0.89788	1.00309	1.09858	1.2044	1.29989	1.4051	1.5012	1.60641	1.7019	1.80711	1.90321	2.00356	2.10451	2.2];
re_7810.pdf = [0.9895	1.21899	1.02616	0.8493	1.00889	0.8493	0.71496	0.80128	0.60186	0.54754	0.33304	0.30817	0.18744	0.13809	0.09809	0.11209	0.06952	0.04147];

plot(re_5210.di_less_dia,re_5210.pdf,'rx',"LineWidth",2,"DisplayName",'Re = 5210 (Exp.)')
hold on
plot(re_7810.di_less_dia,re_7810.pdf,'bx',"LineWidth",2,"DisplayName",'Re = 7810 (Exp.)')
hold on

    fixed_edges = linspace(0.1,3,50);


if (1)%液滴分裂可追踪
    initial_value = 0.2*ones(1,1);
    threshold = data.critical_value(1);
    ratio = [2,18];%固定比例分裂时此比例需小于0.5，随机比例分裂时可以为0-1的任意子区间，gamma随机比例时为输入为[均值,方差],beta分布时为alpha,beta
    tic
    compressed_data = compress_data_process(initial_value, threshold,ratio);
    toc
    dat = reconstruct_data(compressed_data);
    di_less_dat = dat/mean(dat);
    % figure('Visible', 'off'); %创建隐藏图形
    line([threshold/mean(dat) threshold/mean(dat)],[max(re_7810.pdf) min(re_7810.pdf)],'linewidth',2,"DisplayName",'Critical Value');
    h_temp = histogram(di_less_dat,fixed_edges,"Normalization","pdf","DisplayName",'Theory, no Coalesce');
    
    % res.binvalues{1} = h_temp.Values;
    % err = 1;
    % last_bins = 0;
    % iter = 1;
    % while(err>1e-5&&iter<300)
    % 
    %     %开始合并
    %     dat = probabilistic_merge(dat, threshold);
    %     %再次分裂
    %     compressed_data = compress_data_process(dat, threshold,ratio);
    %     dat = reconstruct_data(compressed_data);
    %     di_less_dat = dat/mean(dat);
    %     % line([threshold/mean(dat) threshold/mean(dat)],[max(re_7810.pdf) min(re_7810.pdf)],'linewidth',2);
    %     h1 = histogram(di_less_dat,fixed_edges,"Normalization","pdf","DisplayName",'Theory, with Coalesce');
    % 
    %     err = norm(last_bins-h1.Values);
    %     e(iter) = err;
    %     last_bins = h1.Values;
    %     iter = iter +1;
    %     res.binvalues{iter} = h1.Values;
    %     close(gcf);
    % end
    hold on
    set(gca, 'YScale', 'log');
    legend();
    % ylim([0 10^2]);
    % xlim([0 3]);
    xlabel('$D/\langle D\rangle $','Interpreter','latex');
    ylabel('PDF');
    title(['The number of droplets =',num2str(length(dat))]);
    
    % density_compare(compressed_data, threshold, ratio, 100);%计算PDF并检测收敛
else
    tic
    [dat, iter] = independent_split(1e+6, data.critical_value(1),0.05,0.2);
    toc
    % % 绘制对数直方图
    % figure;
    di_less_dat = dat/mean(dat);
    histogram(di_less_dat,15,"Normalization","pdf","DisplayName",'Theory');
    hold off

    set(gca, 'YScale', 'log');
    legend();
    % ylim([0 10^2]);
    % xlim([0 3]);
    xlabel('$D/\langle D\rangle $','Interpreter','latex');
    ylabel('PDF');
    title(['The number of iter = ',num2str(iter)]);
end


function density_compare(compressed_data,threshold,ratio,max_iter)
    %本函数直接比较尺寸分布的频率
    if size(compressed_data, 2) == 1
        error('输入的必须是压缩数据！')
    end
    load_data = true;%是否读取上次运行的数据


    prev_density = [];
    prev_diameters = [];
    prev_data = [];
    difference = [];
    iter = 1;
    converged = false;
    figure('Position', [100, 100, 1200, 500]);
    
    
    if (1)%从读取数据开始
        load(['data_',num2str(ratio),'\data_split_20.mat'])
        total_num = sum(prev_data(:,2));%粒子总数
        prev_density = prev_data(:,2)/total_num;
        prev_diameters = prev_data(:,1);%所有不同粒子直径
    end

    while ~converged && iter <= max_iter

        total_num = sum(compressed_data(:,2));%粒子总数
        density = compressed_data(:,2)/total_num;
        diameters = compressed_data(:,1);%所有不同粒子直径

        clf;
        if ~isempty(prev_density)
            semilogy(prev_diameters, prev_density, 'bo', 'LineWidth', 1.5, 'DisplayName', '前次频率');
            hold on;
        else
            disp('生成统计区间边界。')
            min_d = min(diameters);
            max_d = max(diameters);
            range = max_d-min_d;
            bin_edge = linspace(min_d-range*0.1,max_d+0.1*range,32);%得到区间边界
            bin_value = diff(bin_edge)/2+bin_edge(1:end-1);%得到区间中点
        end

        
        semilogy(diameters, density, 'rx', 'LineWidth', 1.5, 'DisplayName', '当前频率');% 原始频率图
        % bin_density= bin_counts(diameters, density, bin_edge);%缩减频率图
        % plot(bin_value, bin_density, 'rx', 'LineWidth', 1.5, 'DisplayName', '当前频率');% 缩减后的频率图
        hold off;
        legend('Location', 'best');
        title('频率分布对比');
        xlabel('粒子直径');
        ylabel('频率');
        drawnow;
        fprintf('\n共有%d个不同数据，总数据量为%d\n',length(compressed_data(:,2)),total_num);

        % 计算前后两次频率直方图的差距
        if ~isempty(prev_density)
            difference(iter-1) = norm(prev_density-bin_density,2);
            fprintf('频率分布差值为%.6f',difference(iter-1));
            if difference(iter-1) < 1e-5
                fprintf('频率分布收敛，差值为%.6f',difference(iter-1));
            end
        end


        % 更新前次频率
        prev_density = density;%原始分布
        prev_diameters = diameters;%原始分布
        % prev_density = bin_density;%分区间分布
        % prev_diameters = bin_value;%分区间分布
        prev_data = compressed_data;
        % 再次处理数据

        next_diameters = diameters*10;%倍增
        next_data = [next_diameters,compressed_data(:,2)];%下一步所用的矩阵
        compressed_data = compress_data_process(next_data, threshold,ratio);%下一次分裂

        iter = iter + 1;
        
        filename = fullfile('.',['data_',num2str(ratio)],['data_split_',num2str(iter),'.mat']);
        save(filename,"iter","threshold","compressed_data","difference","prev_data")
    end
end

function [processed_data] = compress_data_process(data, threshold,ratio)
    % 数据处理主函数：将数据按阈值分裂处理，返回（值，计数）矩阵
    % 输入：
    %   data - 原始数据（列向量）
    %   threshold - 分裂阈值
    %   ratio - 分裂比例
    % 输出：
    %   processed_data - 处理后的[Nx2]矩阵，列1为值，列2为对应计数
    
    if(isscalar(ratio))
        split_rule = @(x) deal(x*(1-ratio)^(1/3), x*ratio^(1/3));  %分裂规则
    elseif(length(ratio)==2)
        disp('使用随机分裂。');
    end

    if ~exist(['./data_',num2str(ratio)])
        mkdir(['./data_',num2str(ratio)])
    end

    
    % 初始数据转换
    if size(data,2) == 1
        [values, ~, idx] = unique(data(:)); % 强制转换为列向量
        counts = accumarray(idx, 1);
        current_data = [values, counts];
    else
        current_data = data;
    end

    
    % 迭代处理
    has_changes = true;
    iteration = 0;%进度计算器
    while has_changes
        iteration = iteration + 1;%进度计算器
        % fprintf('正在进行第 %d 次迭代...\n', iteration);

        has_changes = false;
        new_data = [];
        
        for i = 1:size(current_data, 1)
            val = current_data(i, 1);
            count = current_data(i, 2);
            
            if val > threshold
                % ================= 分裂规则区域 ======================
                if(isscalar(ratio))
                    [split_val1, split_val2] = split_rule(val);
                elseif(length(ratio)==2)
                    % 采用均匀分布
                        % ratio_val = ratio(1)+(ratio(2)-ratio(1))*rand(1);
                    % 采用均匀分布

                    % 采用gamma分布,大于1小于0的不考虑
                        % ratio_val = 1;
                        % while (abs(ratio_val-0.5)>=0.5)
                        %     theta_gamma = ratio(2)/ratio(1);
                        %     k_gamma = ratio(1)/theta_gamma;
                        %     ratio_val = gamrnd(k_gamma, theta_gamma, 1, 1);
                        % end
                    % 采用gamma分布,大于1小于0的不考虑

                    % 采用自定义对称beta分布
                        a = ratio(1);
                        b = 5*(val/threshold);%ratio(2);    
                        % y = (betapdf(x, a, b)+betapdf(1-x, a, b))/2;
                        % figure;
                        % plot(x, y, 'LineWidth', 2.5, ...
                        %     'LineStyle', styles{mod(i-1,numel(styles))+1}, ...
                        %     'Color', colors(i,:), ...
                        %     'DisplayName', sprintf('α=%.1f, β=%.1f',a,b));

                        ratio_val = custom_beta_rnd(a,b,1);
                    % 采用自定义对称beta分布

                    split_val1 = val*(1-ratio_val)^(1/3);
                    split_val2 = val*ratio_val^(1/3);
                end
                % ================= 分裂规则结束 ======================
                
                new_data = [new_data; split_val1, count; split_val2, count];
                has_changes = true;
            else
                new_data = [new_data; val, count];
            end

        end
        
        % 合并重复值
        if ~isempty(new_data)
            [merged_values, ~, idx] = unique(new_data(:, 1));
            merged_counts = accumarray(idx, new_data(:, 2));
            current_data = [merged_values, merged_counts];
        else
            current_data = [];
        end

        % 显示当前迭代结果信息
        if isempty(current_data)
            fprintf('迭代 %d 完成，无剩余数据。\n', iteration);
        else
            max_val = max(current_data(:, 1));
            fprintf('迭代 %d 完成，当前数据条目数：%d，最大值：%.5f\n', ...
                    iteration, size(current_data, 1), max_val);
        end
    end
    
    processed_data = current_data;
end

function merged = merge_similar_values(data, tol)
        % 合并相似值（防止值数量膨胀）
        [sorted_val, order] = sort(data(:,1));
        sorted_count = data(order,2);
        
        merge_groups = discretize(sorted_val,...
            [sorted_val(1)-eps; sorted_val(1:end-1)+tol*diff(sorted_val); sorted_val(end)+eps]);
        
        merged_val = accumarray(merge_groups, sorted_val, [], @mean);
        merged_count = accumarray(merge_groups, sorted_count);
        merged = [merged_val, merged_count];
end

function [sumCounts] = bin_counts(values,counts,edges)

    if size(values,2)~=1||size(counts,2)~=1
        error('输入的数据需要将值和数量分离。');
    end

    % 分配数值到区间
    bins = discretize(values, edges);

    
    % 统计有效区间的加权和
    validBins = bins(~isnan(bins));
    validCounts = counts(~isnan(bins));
    numBins = length(edges)-1;
    
    % 使用accumarray处理可能缺失的区间
    sumCounts = accumarray(validBins, validCounts, [numBins,1], @sum, 0);
    
end

function [reconstructed_data] = reconstruct_data(processed_data)

    % 数据重建函数：将（值，计数）格式数据还原为原始数组格式
    % 输入：
    %   processed_data - processDataWithThreshold函数输出的[Nx2]矩阵
    % 输出：
    %   reconstructed_data - 重建后的原始数据格式（列向量）
    
    % 警告提示（适用于大数据时）
    total_count = sum(processed_data(:,2));
    if total_count > 1e5
        warning(['重建数据量达%d，可能占用%.2fMB内存，',...
            '建议仅用于小数据验证，按回车键继续。'], total_count, total_count*8/1e6);
            % pause
    end
    
    % 执行重建
    reconstructed_data = [];
    if ~isempty(processed_data)
        reconstructed_data = repelem(processed_data(:,1), processed_data(:,2));
    end
end

function [dat, iter] = independent_split(n, xu_threshold,p,q)
    %   Inputs:
    %       n            - 液滴个数
    %       p and q      - 上下界
    %       xu_threshold - Threshold value for iterative processing
    %   Outputs:
    %       dat          - Processed data array
    %       iter         - Number of iterations performed

    initial_value = 5;
    % 初始化数组
    dat = initial_value * ones(n, 1);
    iter = 0;

    % 主处理循环
    while sum(dat > xu_threshold) > 1
        % 生成随机数并调整范围到[p, q]
        b = rand(size(dat));
        bb = p + (q-p).*b;
        
        mask = dat > xu_threshold;
        dat(mask) = dat(mask) .* bb(mask);
        
        iter = iter + 1;
    end


end

function new_data = probabilistic_merge(data, threshold)
    %液滴合并
    current_data = data;
    i = 1;

    
    
    while i <= length(current_data)
        % 获取当前元素
        current_element = current_data(i);
        
        % 获取其他元素的索引
        other_indices = setdiff(1:length(current_data), i);
        
        if ~isempty(other_indices)
            % 随机选择一个其他元素
            j = other_indices(randi(length(other_indices)));
            selected_element = current_data(j);
            
            merge_prob = 0.01*((selected_element/threshold)+(current_element/threshold));
            % 根据概率决定是否合并
            if rand() < merge_prob
                % 创建新元素
                merged_element = (current_element^3 + selected_element^3)^(1/3);
                
                % 创建新数组
                mask = true(1, length(current_data));
                mask([i, j]) = false;
                current_data = [current_data(mask); merged_element];
                
                % 重置索引，因为数组结构已改变
                i = max(1, i-1); % 后退一步以处理新元素
            else
                % 不合并则处理下一个元素
                i = i + 1;
            end
        else
            % 没有其他元素可合并时退出
            break;
        end
    end
    
    new_data = current_data;
end







