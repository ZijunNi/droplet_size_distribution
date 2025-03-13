%此程序用于对双对数坐标图像进行处理
%目前的像素位置数据来自于孙超老师实验，下列文献
% Yi, L., Toschi, F., & Sun, C. (2020). Global and local statistics in turbulent emulsions. Journal of Fluid Mechanics, 912.
% 输入已知的像素坐标和对应的绝对坐标
clear,clc;
pixel_coords = [ 630,803;382,803;382,367;]; % N×2矩阵，每行代表一个像素坐标[px, py]
absolute_coords = [ 1e+4,1e-3; 5e+3,1e-3;5e+3,1e-2;]; % N×2矩阵，每行代表对应的绝对坐标[x, y]

% 计算对数绝对坐标
log_abs = log(absolute_coords);

% 拟合x轴转换参数（像素X到log(x)）
X_design = [pixel_coords(:,1), ones(size(pixel_coords,1), 1)];
coeff_x = X_design \ log_abs(:,1); % 最小二乘解

% 拟合y轴转换参数（像素Y到log(y)）
Y_design = [pixel_coords(:,2), ones(size(pixel_coords,1), 1)];
coeff_y = Y_design \ log_abs(:,2);

% 定义转换函数
pixel2abs = @(pixel) deal(...
    exp(coeff_x(1)*pixel(:,1) + coeff_x(2)), ...
    exp(coeff_y(1)*pixel(:,2) + coeff_y(2)));

% 示例使用：将像素坐标[px_input, py_input]转换为绝对坐标
px_input = [399,542,643,787,891,969]; % 输入像素X坐标
py_input = [271,361,431,526,584,621]; % 输入像素Y坐标
% [x_abs, y_abs] = pixel2abs([px_input, py_input]);
% disp(['绝对坐标: (', num2str(x_abs), ', ', num2str(y_abs), ')']);

% 若要批量转换多个坐标，可使用：
test_pixels = [px_input',py_input';];
[x_all, y_all] = pixel2abs(test_pixels);
loglog(x_all,y_all,'x',linewidth=2,DisplayName='Sun Chao Experiment')
hold on
my_data_x = [6.6422e+03,2.5377e+05]';
my_data_y = [1.21e-2,7e-5]';
loglog(my_data_x,my_data_y,'ro',linewidth=2,DisplayName='Our Therory')
hold on
ref = [x_all;my_data_x];
loglog(ref,ref.^(-1.2)*500,'-',DisplayName='$\langle D\rangle\sim Re^{-1.2}$')
hold off
xlabel('$Re$',Interpreter='latex')
ylabel('$\langle D\rangle/d$',Interpreter='latex')
legend(Interpreter="latex")
