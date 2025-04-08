function [z_positions, adjusted_rms, shift_constant, scale_constant] = process_u_rms_profile_corrected(u_velocity_field, z_coordinates)
    % 计算并调整流向速度u的RMS剖面（修正版）
    % 输入:
    %   u_velocity_field - 3D数组 (x,y,z)或(x,z,t)等，包含流向速度分量u
    %   z_coordinates - z方向的坐标数组（空间坐标）
    % 输出:
    %   z_positions - z方向空间坐标
    %   adjusted_rms - 调整后的RMS剖面
    %   shift_constant - 平移常数(使RMS最大值变为0)
    %   scale_constant - 缩放常数(使RMS最小值变为-1)
    
    % 计算原始RMS剖面
    [z_positions, u_rms_profile] = calculate_u_rms_profile(u_velocity_field, z_coordinates);
    
    % 找到原始RMS值的最大值和最小值（注意：是RMS值，不是z坐标）
    max_rms_value = max(u_rms_profile);
    min_rms_value = min(u_rms_profile);
    
    % 计算调整常数
    shift_constant = -max_rms_value;  % 平移量(使RMS最大值变为0)
    scale_constant = -(min_rms_value + shift_constant); % 缩放因子(使平移后的最小值变为-1)
    
    % 应用调整
    shifted_rms = u_rms_profile + shift_constant;
    adjusted_rms = shifted_rms / scale_constant;
    
    % 绘制结果
    figure;
    
    % 子图1: 原始RMS剖面
    subplot(1, 2, 1);
    plot(z_positions, u_rms_profile, 'b-o', 'LineWidth', 2);
    ylabel('u_{rms} (original)');
    xlabel('z position (wall-normal direction)');
    title('Original RMS Profile');
    grid on;
    hold on;
    yline(0, 'k--'); % 添加零线
    hold off;
    
    % 子图2: 调整后的RMS剖面
    subplot(1, 2, 2);
    plot(z_positions, adjusted_rms, 'r-o', 'LineWidth', 2);
    ylabel('u_{rms} (adjusted)');
    xlabel('z position (wall-normal direction)');
    title({'Adjusted RMS Profile', ...
           ['shift=', num2str(shift_constant, '%.4e'), ...
           ', scale=', num2str(scale_constant, '%.4e')]});
    grid on;
    hold on;
    yline(0, 'k--'); % 添加零线
    yline(-1, 'g--'); % 添加-1线
    hold off;
    
    % 统一y轴范围以便比较
    all_rms = [u_rms_profile; adjusted_rms];
    ylim_range = [min(all_rms)-0.1, max(all_rms)+0.1];
    subplot(1, 2, 1); 
    % ylim(ylim_range);
    subplot(1, 2, 2); 
    % ylim(ylim_range);
    
    % 调整图形尺寸和布局
    set(gcf, 'Position', [100, 100, 1200, 500]);
    end
    
