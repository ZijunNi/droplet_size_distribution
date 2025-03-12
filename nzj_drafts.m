%尝试读取数据
clear,clc;
re = 180;
data_fold = ['data_',num2str(re)];
currentFolder = pwd;
data_fold_path = fullfile(currentFolder, data_fold);
data_file_name = ['/Re_tau = ',num2str(re),'.mat'];
load([data_fold_path,data_file_name])
load("self_example_data.mat")

% correction_para = 3;

clc;
close(gcf);
wall_parallel_vel = data.mean_U;%-correction_para;

wall_normal_coors = data.zpos*data.Reynolds_number;
wall_normal_coors = wall_normal_coors(2:end);
wall_parallel_vel = wall_parallel_vel(2:end);

semilogx(wall_normal_coors,wall_parallel_vel,'-x')
hold on
semilogx(wall_normal_coors,log(wall_normal_coors)/0.41+5)
hold on
linear_ref = 1:0.05:32;
semilogx(linear_ref,linear_ref)
hold on
semilogx(linear_ref,log(linear_ref)/0.41+5)
hold on


matching_loc = 1:1:32;
wall_parallel_vel = (U.^2+V.^2);%-correction_para;
wall_parallel_vel = wall_parallel_vel(:,:,2:end);
tic
for i = 1:length(matching_loc)
    res = extract_2d_slice_x_interp(wall_parallel_vel, wall_normal_coors, matching_loc(i)/180,re);
    a(i) = mean(res(:));
end
toc
plot(matching_loc,a,'x',LineWidth=2)
hold off


%%
re_tau = 180;
re = 2.54e5;
delta = 0.005;
r_i = 0.025;
K = 0.74;
(re_tau/(delta*sqrt(K)/r_i))^(1/0.79)


%%
% 输入参数说明：
% data: 三维矩阵，维度为 [a, b, c]
% z: 第三维对应的坐标向量，长度为 c
% z_target: 目标插值坐标（标量）
function result = extract_2d_slice_x_interp(A, z, zq, re_tau)
    [a, b, c] = size(A);
    result = zeros(a, b);
    matching_loc = 12; % 对数区与线性区的匹配位置，应为14.2左右
    save_pos = zq * re_tau; % 转换为y+坐标

    if(max(z)<re_tau)
        z = z*re_tau;% 转换为y+坐标
    end
    
    if(0)%对数插值

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
        result = interpolateLayer(A, z, save_pos, a, b, c);
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



% function result = extract_2d_slice_x_interp(A, z, zq, re_tau)
%     [a, b, c] = size(A);
%     result = zeros(a, b);
%     matching_loc = 0.2; % 对数区与线性区的匹配位置，实际上在14.2左右
%     save_pos = zq * re_tau; % 转换为y+坐标
% 
%     if save_pos > matching_loc
%         % 在对数区进行插值：将z转换为对数坐标
%         z_log = log(z);
%         target = log(save_pos);
%         result = interpolateLayer(A, z_log, target, c);
%     else
%         % 在粘性底层：先插值到匹配位置，再线性插值到目标位置
%         z_log = log(z);
%         target_log = log(matching_loc);
%         result0 = interpolateLayer(A, z_log, target_log, c);
% 
%         % 构造0到matching_loc的线性插值数据
%         z_linear = [0, matching_loc];
%         data_linear = cat(3, zeros(a,b), result0); % 假设y+=0处值为0
%         result = interpolateLayer(data_linear, z_linear, save_pos, 2);
%     end
% end
% 
% function result = interpolateLayer(data, z, target, c)
% 
%     idx = find(z <= target, 1, 'last');
%     if isempty(idx)
%         idx = 1;
%     elseif idx == c
%         idx = c - 1;
%     end
%     % 执行线性插值
%     z1 = z(idx);
%     z2 = z(idx+1);
%     y1 = data(:,:,idx);
%     y2 = data(:,:,idx+1);
%     result = y1 + (target - z1) * (y2 - y1) / (z2 - z1);
% 
% end

