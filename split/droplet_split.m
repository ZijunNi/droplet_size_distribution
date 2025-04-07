% 示例使用
clear,clc
re_tau = 180;
filename = ['../data_',num2str(re_tau),'_old/Re_tau = ',num2str(re_tau),'.mat'];
load(filename)
initial_value = 0.1;
threshold = data.critical_value(1);
ratio = 0.4;%分裂比例小于等于0.5
compressed_data = compress_data_process(initial_value, threshold,ratio);%第一次分裂


density_compare(compressed_data, threshold, ratio, 100);%计算PDF并检测收敛


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
    
    % ratio = 0.5;%小于等于0.5
    split_rule = @(x) deal(x*(1-ratio)^(1/3), x*ratio^(1/3));  %分裂规则
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
    while has_changes
        has_changes = false;
        new_data = [];
        
        for i = 1:size(current_data, 1)
            val = current_data(i, 1);
            count = current_data(i, 2);
            
            if val > threshold
                % ================= 分裂规则区域 ======================
                [split_val1, split_val2] = split_rule(val);
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
            pause
    end
    
    % 执行重建
    reconstructed_data = [];
    if ~isempty(processed_data)
        reconstructed_data = repelem(processed_data(:,1), processed_data(:,2));
    end
end











