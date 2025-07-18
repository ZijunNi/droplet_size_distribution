close all
clear all
clc
load('AEM_data_VEL_FIELDS_progress_save2_Rahul_0.12.mat');

%%
U_ALL = load_data.U_ALL;
V_ALL = load_data.V_ALL;%/utau_est;
W_ALL = load_data.W_ALL;%/utau_est;
Num_layers = 1:6;

n =  2;%累加的层数

% while(n<6)
%     n = n +1;
Re_tau = 100*2^(n-1);

zpos = (0:size(U_ALL,3)-1)*2^(n-1); 
zpos_delta = zpos/(100*2^(n-1));
 
xpos = (0:size(U_ALL,2)-1)*2^(n-1);
xpos_delta = xpos./(100*2^(n-1));

ypos = (0:size(U_ALL,1)-1)*2^(n-1);
ypos_delta = ypos./(100*2^(n-1));

% add up all hier. to get full field
close all

U = sum(U_ALL(:,:,:,1:n),4);
V = sum(V_ALL(:,:,:,1:n),4);
W = sum(W_ALL(:,:,:,1:n),4);


u_mean = mean(U, [1,2], 'omitnan');
w_mean = mean(W, [1,2], 'omitnan');


% 计算脉动量和雷诺应力
for i = 1:length(u_mean)
u_prime(:,:,i) = U(:,:,i) - u_mean(i);
w_prime(:,:,i) = W(:,:,i) - w_mean(i);
end

Reynolds_Shear = mean(u_prime.* w_prime, [1,2]);%, 'omitnan');
profile = squeeze(Reynolds_Shear);
utau_est = sqrt(-min(profile));
% profile = profile/utau_cal^2;
%%
% begin_pos = 2^(7-n)+5;


z_plus = squeeze(z_frame(1,1,:)) * Re_tau;
logical_region = (z_plus > 50) & (z_plus < max(0.1 * Re_tau,60));
indices = find(logical_region);
p_lin = polyfit(log(z_plus(indices)), squeeze(u_mean(indices)), 1);
% p_lin = polyfit(log(squeeze(z_frame(1,1,begin_pos:begin_pos+3))*Re_tau),squeeze(u_mean(begin_pos:begin_pos+3)),1);  % 2:5 corresponds to the range in the mean flow profile I fit too
figure;semilogx(squeeze(z_frame(1,1,:))*Re_tau,squeeze(u_mean),'-x');hold on;
plot([1 Re_tau],polyval(p_lin,log([1 Re_tau])),'--k');
plot([1 Re_tau],polyval([2.5 p_lin(2)],log([1 Re_tau])),'--b');
Uinf_est = 5-p_lin(2);
% get Uinf estimates
Uinf = polyval([2.5 5],log(Re_tau));
% U_ALL = U_ALL +Uinf;

%%
% add up all hier. to get full field
close all

U = sum(U_ALL(:,:,:,1:n),4)/utau_est+Uinf;
V = sum(V_ALL(:,:,:,1:n),4)/utau_est;
W = sum(W_ALL(:,:,:,1:n),4)/utau_est;


u_mean = mean(U, [1,2], 'omitnan');
w_mean = mean(W, [1,2], 'omitnan');


%% 绘制雷诺应力对比图
close all

figure;
sgtitle(['$Re_\tau = $',num2str(Re_tau)],'Interpreter', 'latex');
% 设置全局边距（单位为归一化坐标，范围[0,1]）
% [左，下，右，上]的边距
bottom_margin = 0.15; % 增加底部边距
subplot_margin = 0.05; % 子图之间的间距

% 调整子图位置
subplot(2,2,1);
set(gca, 'Position', [0.1, 0.7, 0.35, 0.2]); % [左，下，宽，高]
semilogx(zpos_delta*Re_tau, squeeze(u_mean));
hold on 
semilogx(zpos_delta*Re_tau,log(zpos_delta*Re_tau)/0.41+5,DisplayName='Log Law');
xlabel('Wall-normal Distance $z^+$', 'Interpreter', 'latex');
ylabel('Raw $\overline{u}$', 'Interpreter', 'latex');

subplot(2,2,2);
set(gca, 'Position', [0.55, 0.7, 0.35, 0.2]);
semilogx(zpos_delta*Re_tau, squeeze(w_mean));
xlabel('Wall-normal Distance $z^+$', 'Interpreter', 'latex');
ylabel('Raw $\overline{w}$', 'Interpreter', 'latex');

subplot(2,2,3);
z = zpos_delta;
set(gca, 'Position', [0.1, 0.2, 0.35, 0.4]);
plot(z*Re_tau, profile, 'LineWidth', 1.5);
ylim([1.1*min(profile), 0]);
xlabel('Wall-normal Distance $z^+$', 'Interpreter', 'latex');
ylabel('$\langle u''w'' \rangle$', 'Interpreter', 'latex');
grid on;

subplot(2,2,4);
x = zpos_delta; 
y = squeeze(u_mean);
% 计算中点处的导数
[dy_dx, x_mid] = central_diff_midpoints(y, x);

set(gca, 'Position', [0.55, 0.2, 0.35, 0.4]);
plot(x_mid*Re_tau, dy_dx);
xlabel('Wall-normal Distance $z^+$', 'Interpreter', 'latex');
ylabel(' $\frac{d\overline{u}}{dz}$', 'Interpreter', 'latex');
grid on;

% 添加底部文字（使用 annotation）
text_content = {
    ['Minimum value of Reynolds stress profile = ', num2str(min(profile))], 
    ['Estimated value of u_tau = ', num2str(sqrt(-min(profile)))]};

annotation('textbox', [0.1, 0.01, 0.8, 0.1], ... % [左，下，宽，高]
    'String', text_content, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 10, ...
    'Interpreter', 'none');


%% 将平均场替换为满足law-of-wall的场

U_prime = (U - u_mean); % 脉动分量 [Nx×Ny×Nz]

u_mean_new = linspace(0,Re_tau,101);
u_mean_new = reshape(u_mean_new, 1, 1, []);
U_new = u_mean_new + U_prime; % [Nx×Ny×Nz]
figure;imagesc(xpos_delta,ypos_delta,squeeze(U_new(:,:,2))'); daspect([1 1 1]); %  wall parallel plane 
a = squeeze(U_new(:,:,2))';
mean(a(:));
%%
U = U_new;
save(['converted_de_silva_data_re_tau_',num2str(Re_tau),'.mat'],"U","V","W","xpos_delta","ypos_delta","zpos_delta");

%% 绘制处理后的速度剖面对比图
% close all
% 
% figure;
% sgtitle(['$Re_\tau = $',num2str(Re_tau)],'Interpreter', 'latex');
% bottom_margin = 0.15; % 增加底部边距
% subplot_margin = 0.05; % 子图之间的间距
% 
% semilogx(zpos_delta*Re_tau, squeeze(u_mean),'-x',DisplayName=['Velocity profile, $Re_\tau = $',num2str(Re_tau)],LineWidth=2);
% hold on 
% semilogx(zpos_delta*Re_tau,log(zpos_delta*Re_tau)/0.41+5,DisplayName='Log Law, $u^+=2.5\log(y^+)+5$');
% hold off
% 
% % 创建文本标注内容（使用sprintf格式化输出）
% text_content = {
%     sprintf('Utau =   %.6f', utau_est),  % 第一行
%     sprintf('Uinf  = %.6f', Uinf)   % 第二行
% };
% 
% % 在左上角添加文本标注（使用归一化坐标定位）
% text(0.03, 0.95, text_content,...
%     'Units', 'normalized',...         % 使用归一化坐标系
%     'VerticalAlignment', 'top',...    % 顶部对齐
%     'HorizontalAlignment', 'left',... % 左对齐
%     'BackgroundColor', [1 1 1],...    % 白色背景
%     'EdgeColor', 'k',...              % 黑色边框
%     'Margin', 3,...                   % 文本边距
%     'FontSize', 10,...                % 字号大小
%     'FontName', 'Consolas');          % 等宽字体更美观
% 
% % 调整坐标区域使文本不被遮挡（可选）
% set(gca, 'Position', [0.15 0.15 0.75 0.75]);
% 
% legend(Interpreter="latex",Location="southeast")
% xlabel('Wall-normal Distance $z^+$', 'Interpreter', 'latex');
% ylabel('$\overline{u^+}$', 'Interpreter', 'latex');
% 
% % 
% % subplot(1,2,2);
% % set(gca, 'Position', [0.6, 0.1, 0.4, 0.8]);
% % semilogx(zpos_delta*Re_tau, squeeze(w_mean));
% % xlabel('Wall-normal Distance $z^+$', 'Interpreter', 'latex');
% % ylabel('$\overline{w^+}$', 'Interpreter', 'latex');
% 
% 
% %% 保存图像
%     figure_name = ['Re_tau = ',num2str(Re_tau),'.pdf'];
%     filename = fullfile(pwd,figure_name);
%     exportgraphics(gcf, filename, 'ContentType', 'vector');
%     % close all;
