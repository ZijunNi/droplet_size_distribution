function [threshold_series,time_m] = threshold_line(u,v,dx,targetProbability,height_distance,zpos,Re)

    % 计算wall-normal速度平方场
    wall_normal_velocity = (u.^2+v.^2);
    interped_2D_field = extract_2d_slice_x_interp(wall_normal_velocity,zpos,height_distance,Re);
    % 计算开始和结束的区间位置
    full_data = sort(reshape(interped_2D_field,[],1));%一维的完整数据
    [min_threshold,max_threshold] = calculate_percentiles(full_data,0.001);
    
    numSignals = length(interped_2D_field(1,:)); % 信号组数，与展向采样点数量相同
    sampling_interval = dx; % 采样间隔，与流向采样点（时间或空间）间距相同
    
    num_threshold = 100;% 阈值采样点数量
    threshold_series = linspace(min_threshold,max_threshold,num_threshold);%所需要计算的阈值序列
    
    % 开始计算
    % ii = 1;
    time_m = zeros(1,num_threshold);

    parfor ii=1:num_threshold%此处可引入并行加快速度

        threshold = threshold_series(ii);
        % 用于存储超过阈值的序列长度
        % sequenceLengths = zeros(numSignals, 1); % 假设最坏情况长度为 numSignals
        sequenceLengths = [];
        % 遍历每组信号
        for i = 1:numSignals
            signal = interped_2D_field(:,i); % 获取当前信号

            [duration,~] = calculate_time_above_threshold(signal,threshold,sampling_interval);

            sequenceLengths = [sequenceLengths;duration];
        end
        [~,time_m(ii)] = calculate_percentiles(sequenceLengths,targetProbability);%真实统计数据
    
    end
end