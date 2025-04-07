function [z_positions, u_rms_profile] = calculate_u_rms_profile(u_velocity_field, z_coordinates)
    % 计算流向速度u沿z方向的RMS剖面
    [~, ~, nz] = size(u_velocity_field);
    u_rms_profile = zeros(nz, 1);
    z_positions = z_coordinates;
    
    for k = 1:nz
        u_plane = u_velocity_field(:, :, k);
        u_mean = mean(u_plane(:), 'omitnan');
        u_fluc = u_plane - u_mean;
        u_rms = sqrt(mean(u_fluc(:).^2, 'omitnan'));
        u_rms_profile(k) = u_rms;
    end
    end