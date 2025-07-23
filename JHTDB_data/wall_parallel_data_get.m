% 提取二维平面（平行于壁面）数据
 
clear all;
close all;
clear,clc
% ---- Enter user JHTDB token ----
authkey = 'cn.edu.pku.stu.nizijun-955f3869';  
% the above is a default testing token that works for queries up to 4096 points
% for larger queries, please request token at Please send an e-mail to 
% turbulence@lists.johnshopkins.edu including your name, email address, 
% and institutional affiliation and department, together with a short 
% description of your intended use of the database.
%
% ---- select dataset ----
dataset =  'channel';


% ----- Initialize getData parameters (except time and points) -----
variable = 'velocity';
temporal_method = 'none'; 
spatial_method = 'lag4';
spatial_operator  = 'field';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D plane demo points : evenly spaced over a 2D plane lying along one of the primary axes
%     - time : the time to be queried (snapshot number for datasets without a full time evolution).
%     - nx, nz : number of points along each axis. total number of points queried will be n_points = nx * nz.
%     - x_points, y_points, z_points : point distributions along each axis, evenly spaced over the specified ranges.
%     - linspace(axis minimum, axis maximum, number of points).
%     - points : the points array evenly spaced out over the 2D plane.
%     - points array is instantiated as an empty array that will be filled inside the for loops.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_pos = linspace(0.1,3,10)/1000;% 提取y^+在0.1到3之间的10个平面

 for num_pos=1:10
    tic
    time = 10;
    
    nx = 2048;
    nz = 400;
    n_points = nx * nz;
     
    points = zeros(n_points,3);
    
    x_points = linspace(0.0, 8 * pi, nx);
    y_points = 1-y_pos(num_pos);
    z_points = linspace(0.0, 0.8 * pi, nz);
    
    for i = 1 : nx
        for j = 1 :nz
           points(j +(i - 1) * nz, 1) = x_points(i);
           points(j +(i - 1) * nz, 2) = y_points;
           points(j +(i - 1) * nz, 3) = z_points(j);
        end
    end
    
    % ---- GetData ----
    fprintf('\nRequesting %s at %i points...\n', variable, n_points);
    result = getData(authkey, dataset, variable, time, temporal_method, spatial_method, spatial_operator, points);
    
    % 以下为绘图部分
    if (nx >= 2) & (nz >= 2)
        % which component (column) of the data to plot (1-based index, so the first component is specified as 1).
        plot_component = 1;
        % 
        % % ---- Display sample results on screen ----
        % figure1 = figure('Color', [1 1 1], 'InvertHardcopy', 'off', 'PaperSize', [20.98 29.68]);
        % axes1 = axes('FontSize', 16, 'LineWidth', 1.5, 'Parent', figure1, ...
        %     'XScale', 'lin', 'YScale', 'lin', 'Position', [0.18 0.18 0.76 0.76]); 
        % box(axes1, 'on'); 
        % hold(axes1, 'all'); 
        % 
        % % Plotting data
        results = reshape(result(:,plot_component), [nz, nx]); 
        % contourf(axes1, x_points, z_points, results, 300, 'LineColor','none');
        % set(axes1, 'YDir', 'normal');
        % colormap('hot')
        % 
        % % Title and labels
        % title([dataset, ' ', '(', variable, spatial_operator, ')'], 'FontSize', 15);
        % xlabel(axes1, 'X');
        % ylabel(axes1, 'Z');
        % colorbar('FontSize', 16, 'Parent', figure1);
        % set(axes1, 'DataAspectRatio', [1 1 1]); 
        % axis tight; 
        % set(axes1, 'XTickLabel', axes1.XTick, 'YTickLabel', axes1.YTick); 
        % set(axes1, 'TickDir', 'out', 'TickLength', [0.02 0.02]);
        pause(1);
    end
    
    disp(mean(results(:)));
    save(['data_JHTDB_y_plus_',num2str(y_pos(num_pos)*1000),'.mat'],"result","x_points","y_points",'z_points');
    toc
 end
     %提示音
    load chirp 
    sound(y,Fs)