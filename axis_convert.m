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


% DNS datasets 
my_data_x = [6.6422e+03,58211.0249469943]';
my_data_y = [1.32766e-2,10.9e-4]';
% de Silva Datasets
de_silva_data_x = [7589.84818005749	18250.8760745211	43886.8432655531	105532.194945100]';
de_silva_data_y = [1.08652e-2 2.13461e-3 4.36815e-4 9.36184e-5]';


ref = [x_all;58211.0249469943];
loglog(ref,ref.^(-1.2)*500,'--','color', [0.5 0.5 0.5],DisplayName='$\langle D\rangle\sim Re^{-1.2}$',LineWidth=1);hold on
loglog(de_silva_data_x,de_silva_data_x.^(-1.84)/de_silva_data_x(1)^(-1.84)*de_silva_data_y(1),'-.','color', [0.5 0.5 0.5], ...
    DisplayName='$\langle D\rangle\sim Re^{-1.84}$',LineWidth=1);

loglog(x_all,y_all,'kx',linewidth=2,DisplayName='Experiment from Yi $et\, al.$ 2022 ',MarkerSize=10)

loglog(my_data_x,my_data_y,'ks',linewidth=2,DisplayName='DNS Dataset',MarkerSize=10)
loglog(de_silva_data_x,de_silva_data_y,'kd',linewidth=2,DisplayName='AEH fields from de Silva $et\, al.$ 2016',MarkerSize=10);
xlim([4000 70000])
ylim([1e-4 0.05])
hold off
xlabel('$Re$',Interpreter='latex')
ylabel('$\langle D\rangle/d$',Interpreter='latex')
legend(Interpreter="latex")
daspect([1/(0.05-1e-4) 1/(70000-4000) 1])