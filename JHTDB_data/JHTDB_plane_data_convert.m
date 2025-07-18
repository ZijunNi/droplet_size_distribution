% 读取多个平面速度场
% results = reshape(result(:,plot_component), [nz, nx]); % 将原始数据转化为单分量矩阵

clear,clc
u_tau = 0.0499;
Re_tau = 1000;
%% 读取部分
    % 定义特定字符串和文件夹路径
    specific_str = 'data_JHTDB_y_plus_'; % 请替换为你的特定字符串
    folder_path = './'; % 请替换为你的文件夹路径
    
    % 转义特定字符串中的正则表达式特殊字符
    specific_str_escaped = regexptranslate('escape', specific_str);
    
    % 构造匹配文件名的正则表达式模式
    pattern = ['^', specific_str_escaped, '(\d+\.\d*|\.\d+)\.mat$'];
    
    % 获取文件夹内所有.mat文件
    file_list = dir(fullfile(folder_path, '*.mat'));
    
    % 初始化结构数组存储匹配文件信息
    matched_files = struct('name', {}, 'number', {});
    
    % 遍历文件并筛选符合要求的文件
    for i = 1:numel(file_list)
        file_name = file_list(i).name;
        tokens = regexp(file_name, pattern, 'tokens');
        
        if ~isempty(tokens)
            % 提取浮点数字符串并转换为数值
            num_str = tokens{1}{1};
            num_val = str2double(num_str);
            
            % 验证转换是否有效
            if isnan(num_val)
                warning('无法转换数值字符串: %s', num_str);
                continue;
            end
            
            % 存储文件信息
            matched_files(end+1) = struct(...
                'name', fullfile(folder_path, file_name), ...
                'number', num_val ...
            );
        end
    end
    
    % 按浮点数值对文件排序
    [~, sort_idx] = sort([matched_files.number]);
    matched_files = matched_files(sort_idx);
    
    % 加载所有匹配文件的数据
    for i = 1:numel(matched_files)
        file_data = load(matched_files(i).name);
        % 在此处处理数据
        nz = length(file_data.z_points);
        nx = length(file_data.x_points);
        file_data.U = reshape(file_data.result(:,1), [nz, nx]);% 流向
        file_data.V = reshape(file_data.result(:,3), [nz, nx]);% 垂向
        file_data.W = reshape(file_data.result(:,2), [nz, nx]);% 展向
        file_data.y_points = 1-file_data.y_points;
        file_data = rmfield(file_data, 'result');
        % 例如: 将数据存储到工作区或结构数组中
        fprintf('已加载文件: %s (数值 = %.2f)\n', matched_files(i).name, matched_files(i).number);
        % 示例: 将每个文件的数据存入元胞数组
        all_data{i} = file_data; 
    end
% 读取结束

%% 速度剖面验证
for i = 1:length(all_data)

    y_pos(i) = all_data{i}.y_points*Re_tau;

    mean_U(i) = mean(all_data{i}.U(:))/u_tau;

end
figure;
plot(y_pos,mean_U,'-x',DisplayName='Velocity Profile');
hold on
plot(y_pos,y_pos,DisplayName='Linear Law')
hold off
%% 速度场合并

% 初始化
U = [];
V = [];
W = [];
y_points = [];

% 开始合并
for i = 1:length(all_data)
    U = cat(3,U,all_data{i}.U);
    V = cat(3,V,all_data{i}.V);
    W = cat(3,W,all_data{i}.W);
    y_points = [y_points all_data{i}.y_points];
end
% 转化为无量纲速度
U = U/u_tau;
V = V/u_tau;
W = W/u_tau;

x_points = all_data{i}.x_points;
z_points = all_data{i}.z_points;

%% 速度场合并验证
clear "mean_U"

figure;mean_U = squeeze(mean(mean(U,1),2));plot(y_points*Re_tau,mean_U);% velocity profile

figure;imagesc(x_points,z_points,squeeze(U(:,:,10))'); %daspect([1 1 1]); %  wall parallel plane 


%% 速度场格式转化

% post_processing.m 数据格式如下：
% i = Streamwise, j = Spanwise, k = Wall-normal 
% U = Streamwise Velocity
% V = Spanwise Velocity
% W = Wall-normal Velocity

% file_data.U = reshape(file_data.result(:,1), [nz, nx]);% 流向
% file_data.V = reshape(file_data.result(:,3), [nz, nx]);% 垂向
% file_data.W = reshape(file_data.result(:,2), [nz, nx]);% 展向

% 速度场方向转化
U = U;%流向速度仍为U
mid = V;
V = W;%V为展向速度
W = mid;%W为垂向速度

% 速度场坐标轴方向转化
U = permute(U,[2,1,3]);
V = permute(V,[2,1,3]);
W = permute(W,[2,1,3]);
% 坐标数据方向转化
xpos_delta = x_points;
ypos_delta = z_points;
zpos_delta = y_points;

%% 保存数据为post_processing格式

save("../data/data_JHTDB.mat","U","V","W","xpos_delta","ypos_delta","zpos_delta");