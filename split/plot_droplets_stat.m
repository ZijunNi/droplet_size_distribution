function [edges,counts,errors] = plot_droplets_stat(history,M,S,plot_switch)
% 本函数处理coalescence_breakup.m得到的液滴迭代序列
% history - 液滴迭代序列
% M - 统计液滴样本数
% S - 统计液滴样本间距
% plot_switch - 完整统计信息绘制开关，1为打开

% errors - PDF与实验的残差，默认是0
    errors = 0;


    % 计算必要统计量
    num_droplets = cellfun(@length, history);     % 每个时间步的液滴数量
    history_mean_size = cellfun(@mean, history);  % 每个时间步的液滴平均尺寸
    history_mean_std = cellfun(@std, history);    % 每个时间步的液滴尺寸标准差
    
    % 提取最后200步（每隔5步）的液滴尺寸数据
    N = length(history);
    start_idx = max(1, N - S*(M-1));  % 防止索引小于1
    selected_history = history(start_idx:S:N);
    final_sizes = cell2mat(selected_history');
    

    % 自定义分箱边界
    % edges = [0.02:0.12:1.04, 1.05:0.1:3, inf];%绘图边界
    edges = [0.02:0.02:3, inf];%检查边界
    

    if(plot_switch)%绘制完整统计信息图
        % 创建主图形窗口
        % figure;
        
        % ===== 左侧大图：尺寸分布 =====
        subplot(2, 2, [1, 3]);
        normalized_final_sizes = final_sizes / mean(final_sizes);
            
        % 计算并绘制PDF
        [counts, edges] = histcounts(normalized_final_sizes, edges, "Normalization", "pdf");
        centers = edges(1:end-1) + diff(edges)/2;
        plot(centers, counts, '-gx', 'LineWidth', 1.5, 'DisplayName', 'Simulation');

        hold on;
    
        % ===== 尾部对数线性拟合 =====
        % 识别尾部数据范围（从最大值到结束）
        [~, maxIdx] = max(counts);
        endIdx = find(centers>2.2,1);
        tailIdx = maxIdx:endIdx;
        
        % 提取尾部数据（排除零或负值）
        x_fit = centers(tailIdx);
        y_fit = counts(tailIdx);
        validIdx = y_fit > 0;  % 仅使用正值
        x_fit = x_fit(validIdx);
        y_fit = y_fit(validIdx);
        
        if length(x_fit) > 1  % 确保有足够的数据点
            % 进行对数线性拟合（对y取对数后拟合）
            logy = log(y_fit);
            p = polyfit(x_fit, logy, 1);
            slope = p(1);
            intercept = p(2);
            
            % 生成拟合曲线数据
            x_fit_line = linspace(min(x_fit), max(x_fit), 100);
            y_fit_line = exp(polyval(p, x_fit_line));
            
            % 绘制拟合曲线
            plot(x_fit_line, y_fit_line, 'r--', 'LineWidth', 2, 'DisplayName', 'Exponential Fit');
            
            % 添加斜率信息到图例
            slope_str = sprintf('Slope: %.4f', slope);
            text(0.6, 0.7, slope_str, 'Units', 'normalized', 'FontSize', 12, ...
                 'BackgroundColor', [1 1 1 0.7]);
            
        else
            warning('尾部数据不足或无效，跳过拟合');
        end

        
        % 实验数据
        if (1)
            re_5210.di_less_dia = [0.19847	0.39978	0.59622	0.80178	1.00309	1.2044	1.4051	1.60641	1.80286	2.00842	2.2	2.40131	2.60201	2.80332];
            re_5210.pdf = [0.03176	0.43487	1.08767	1.00889	0.72896	0.51658	0.48737	0.24535	0.18075	0.12084	0.05521	0.03559	0.0234	0.02726];
            
            re_7810.di_less_dia = [0.50013	0.60108	0.7063	0.80178	0.89788	1.00309	1.09858	1.2044	1.29989	1.4051	1.5012	1.60641	1.7019	1.80711	1.90321	2.00356	2.10451	2.2];
            re_7810.pdf = [0.9895	1.21899	1.02616	0.8493	1.00889	0.8493	0.71496	0.80128	0.60186	0.54754	0.33304	0.30817	0.18744	0.13809	0.09809	0.11209	0.06952	0.04147];
            
            
            plot(re_5210.di_less_dia, re_5210.pdf, 'rs', "LineWidth", 2, "DisplayName", 'Re = 5210 (Exp.)');
            plot(re_7810.di_less_dia, re_7810.pdf, 'bx', "LineWidth", 2, "DisplayName", 'Re = 7810 (Exp.)');
            
            % ===== 计算残差 ===== 
            
            x_target = [re_5210.di_less_dia re_7810.di_less_dia];
            y_target = [re_5210.pdf re_7810.pdf];
            
            
            [errors, in_range] = calcLogResiduals(centers, counts, x_target, y_target);
            
            str1 = sprintf('最小绝对残差: %.4f log-units\n', min(errors(in_range)));
            str2 = sprintf('最大绝对残差: %.4f log-units\n', max(errors(in_range)));
            str3 = sprintf('平均绝对残差: %.4f log-units\n', mean(errors(in_range)));
            text(0.1, 0.1, {str1,str2,str3}, 'Units', 'normalized', 'FontSize', 12, ...
            'BackgroundColor', [1 1 1 0.7]);
        end
        hold off;
        
        % 设置坐标轴
        set(gca, 'YScale', 'log');
        legend('Location', 'northeast');
        title(sprintf('Final Droplet Size Distribution \n Mean Size = %.5e and %6d Droplets', mean(final_sizes),length(final_sizes)));
        xlabel('Normalized Droplet Size $D/\langle D\rangle$', 'Interpreter', 'latex');
        ylabel('Probability Density');
        
        % ===== 右上图：液滴数量变化 =====
        subplot(2, 2, 2);
        plot(1:length(num_droplets), num_droplets, '-o');
        xlabel('Time Steps');
        ylabel('Number of Droplets');
        grid on;
        title('Droplet Count Evolution');
        
        % ===== 右下图：统计量变化 =====
        subplot(2, 2, 4);
        semilogy(1:length(num_droplets), history_mean_size, '-x', 'DisplayName', 'Mean Size');
        hold on;
        semilogy(1:length(num_droplets), history_mean_std, '-o', 'DisplayName', 'Size Std. Dev.');
        hold off;
        
        % 设置图形属性
        legend('Location', 'best');
        title('Statistical Metrics');
        xlabel('Time Steps');
        ylabel('Metric Value');
        grid on;
        
        % 调整布局
        set(gcf, 'Position', [100, 100, 900, 600]);
    else%仅绘制PDF
        normalized_final_sizes = final_sizes / mean(final_sizes);
        
        
        % 计算并绘制PDF
        [counts, edges] = histcounts(normalized_final_sizes, edges, "Normalization", "pdf");
        centers = edges(1:end-1) + diff(edges)/2;
        plot(centers, counts,'-x', 'LineWidth', 1.5);
    end

    if(0)% 保存统计量
        stat_fold = 'stat';
        data_of_droplets = final_sizes;
        date_str = datestr(now);
        % 需要替换的字符列表
        chars_to_replace = [':', ' ','-'];
        replace_positions = ismember(date_str, chars_to_replace);
        % 执行替换
        date_str(replace_positions) = '_';
        save(fullfile(stat_fold,[date_str,'.mat']),"counts","centers","data_of_droplets")
    end
end