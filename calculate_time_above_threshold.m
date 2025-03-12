function [durations,crossings] = calculate_time_above_threshold(signal, threshold, sampling_interval)
    % 计算周期性信号中超过特定阈值的时间长度
    %
    % 输入参数:
    % signal: 信号序列 (向量)
    % threshold: 阈值 (标量)
    % sampling_interval: 采样间隔 (标量)
    %
    % 输出参数:
    % time_above_threshold: 超过阈值的时间长度 (标量)

    num_samples = length(signal); % 信号长度
    durations = []; % 存储每个时间段的长度
    crossings = []; % 存储每个时间段的开始和结束时间

    % 检查阈值是否小于信号的最小值或大于最大值
    if threshold > max(signal)
        % 如果阈值大于最大值，信号始终超过阈值
        durations = [];
        crossings = [];
        return;
    elseif threshold < min(signal)
        % 如果阈值小于最小值，信号从未超过阈值
        
        durations = num_samples * sampling_interval; % 整个周期都超过阈值
        crossings = [0, num_samples * sampling_interval]; % 整个周期的时间段
        return;
    end


    is_above_threshold = false; % 标记当前是否超过阈值
    start_time = 0; % 初始化开始时间

    

    % 遍历信号序列，考虑周期性
    for i = 1:num_samples
        % 当前点和下一个点（考虑周期性）
        current_sample = signal(i);
        next_sample = signal(mod(i, num_samples) + 1); % 使用 mod 处理周期性

        % 检查是否从下往上跨越阈值（进入超过阈值的区间）
        if(current_sample <= threshold && next_sample >= threshold)
            % 使用线性插值计算精确的跨越时间
            t_cross = (threshold - current_sample) / (next_sample - current_sample) * sampling_interval + (i-1) * sampling_interval;
            start_time = t_cross; % 记录开始时间
            is_above_threshold = true; % 标记为超过阈值
        end

        % 检查是否从上往下跨越阈值（离开超过阈值的区间）
        if(current_sample >= threshold && next_sample <= threshold)
            % 使用线性插值计算精确的跨越时间
            t_cross = (threshold - current_sample) / (next_sample - current_sample) * sampling_interval + (i-1) * sampling_interval;
            end_time = t_cross; % 记录结束时间
            is_above_threshold = false; % 标记为未超过阈值

            % 计算时间段的长度并保存
            duration = end_time - start_time;
            if(start_time>end_time)
                disp('ERROR! Duration is less than ZERO!')
            end
            durations = [durations; duration];
            crossings = [crossings; start_time, end_time];
        end
    end

    % 如果信号最后一个点仍然超过阈值，处理周期结束的情况
    if is_above_threshold && signal(1) >=threshold
        end_time = num_samples * sampling_interval + crossings(1,2); % 周期结束时间+第一个序列结束的时间
        duration = end_time - start_time;
        durations = durations(2:end, :);%删去第一行
        crossings = crossings(2:end, :);%删去第一行
        durations = [durations; duration];
        crossings = [crossings; start_time, end_time];
    end
end