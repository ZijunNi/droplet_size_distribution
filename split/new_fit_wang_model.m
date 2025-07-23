clear,clc
i=2;

filename = ['auto_saved_data_d_',num2str(0.15*i),'_eps_24.76.mat'];

load(filename)

xdata = result(:,1);
ydata = result(:,2);

% 假设已有数据 xdata 和 ydata
ft = fittype(@(a,b,c,k,x) a*(1-exp(-b*x)).*(c*x.^2+k*x), ...
    'independent', 'x');
% 初始参数估算（根据数据特征调整）
y_max = max(ydata);
a0 = y_max / 0.8;     % 高度系数
b0 = 5 / 0.1;         % 上升段系数
c0 = 10 / (0.05)^2;   % 峰宽系数 
k0 = 0.02;             % 下降斜率

fo = fitoptions('Method', 'NonlinearLeastSquares', ...
               'StartPoint', [a0, b0, c0, k0]);

    %                'Lower', [0, 20, 100, -20],...  % 参数下界
    % 'Upper', [20, 200, 1000, ]);   

[fitresult, gof] = fit(xdata, ydata, ft, fo);
plot(fitresult, xdata, ydata);