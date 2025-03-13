function result = extract_2d_slice_x_interp(A, z, zq, re_tau)
%zq为物理空间中的插值位置，z为物理空间中的垂向坐标点位置，都不是无量纲量
    [a, b, c] = size(A);
    result = zeros(a, b);
    matching_loc = 12; % 对数区与线性区的匹配位置，应为14.2左右
    save_pos = zq * re_tau; % 转换为y+坐标

    % if(max(z)<re_tau)
    %     z = z*re_tau;% 转换为y+坐标
    % end
    
    if(1)%对数插值

        if(save_pos > matching_loc)
            % 在对数区进行插值：将z转换为对数坐标
            z_log = log(z);
            target = log(save_pos);
            result = interpolateLayer(A, z_log, target, a, b, c);
        else
            % 在粘性底层：先插值到匹配位置，再线性插值到目标位置
            z_log = log(z);
            target_log = log(matching_loc);
            result0 = interpolateLayer(A, z_log, target_log, a, b, c);
    
            % 构造0到matching_loc的线性插值数据
            z_linear = [0, matching_loc];
            data_linear = cat(3, zeros(a,b), result0); % 假设y+=0处值为0
            result = interpolateLayer(data_linear, z_linear, save_pos, a, b, 2);
        end

    else%完全线性插值
        z = [0;z];
        A = cat(3, zeros(a,b), A);
        result = interpolateLayer(A, z, zq, a, b, c);
    end
end

function result = interpolateLayer(data, z, target, a, b, c)

    idx = find(z <= target, 1, 'last');
    if isempty(idx)
        idx = 1;
    elseif idx == c
        idx = c - 1;
    end
    % 执行线性插值
    z1 = z(idx);
    z2 = z(idx+1);
    y1 = data(:,:,idx);
    y2 = data(:,:,idx+1);
    result = y1 + (target - z1) * (y2 - y1) / (z2 - z1);

end