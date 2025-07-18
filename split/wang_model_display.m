clear,clc
size = [3,6];
epsilon = [1,2,4];
% 假设数据存储结构：X 为横坐标，Y 为 6 组纵坐标数据（每组长度相同）
% 参数 A 有 m 个取值，参数 B 有 n 个取值，总组数 = m × n = 6
m = 2; % A 参数取值个数（示例）
n = 1; % B 参数取值个数（示例）

% 1. 定义颜色方案：为每个A值分配一个基色，同A值下用深浅区分B值
baseColors = [0.0, 0.3, 0.9;  % A1: 蓝色系 (RGB)
              0.9, 0.2, 0.1];  % A2: 红色系 (RGB)

% 2. 定义线型方案：每个B值分配唯一线型
lineStyles = {':', '--', '-.'}; % 线型库，确保数量 >= n

% 3. 生成颜色深浅梯度（同A值下不同B的深浅）
shadeFactors = linspace(0.5, 1, n); % 从深到浅（0.5暗 -> 1亮）

% 4. 绘图
figure; hold on;
legendEntries = cell(m*n, 1); % 存储图例文本

for iA = 1:m       % 遍历A参数
    baseColor = baseColors(iA, :); % 当前A的基色
    
    for iB = 1:n   % 遍历B参数
        % 计算当前曲线索引 (按A优先排列)
        idx = (iA-1)*n + iB; 
        
        % 生成深浅颜色：基色 × 明暗系数
        curveColor = baseColor * 1;%shadeFactors(iB);
        
        % 选择线型 (循环使用线型库)
        lineStyle = lineStyles{mod(iA-1, numel(lineStyles)) + 1};

        filename = ['ddsd_wang_',num2str(size(iA)),'_',num2str(epsilon(iB)),'.mat'];

        load(filename)

        % result = smoothCurve(result);
        result = [0,0;result(1:20,:); smoothCurve(result(25:end,:)); result(end,:)];
        
        % 绘图
        plot(result(:,1),result(:,2), ...
            'Color', 'black', ...
            'LineStyle', lineStyle, ...
            'LineWidth', 1, ...
            DisplayName=['$\varepsilon = ',num2str(size(iA)),'\mathrm{m^3/s^2}, d = ',num2str(epsilon(iB)),'\rm{mm}$'] );
        
        % 记录图例标签（示例格式）
        % legendEntries{idx} = sprintf('A%d, B%d', iA, iB);
    end
end
legend(Interpreter='latex');
xlabel('$f_v$',Interpreter='latex')
ylabel('$\beta(f_v;d,\varepsilon)$',Interpreter='latex')
xlim([0 0.5]);ylim([0 10])
daspect([1/10 1/0.5 1])
grid on
%%

for i=2:4
        % 选择线型 (循环使用线型库)
        lineStyle = lineStyles{mod(i, numel(lineStyles)) + 1};

        filename = ['auto_saved_data_d_',num2str(0.15*i),'_eps_24.76.mat'];

        load(filename)
        % result = smoothCurve(result(1:end,:));
        % 绘图
        plot(result(:,1),result(:,2), ...
            'Color', 'black', ...
            'LineStyle', lineStyle, ...
            'LineWidth', 1, ...
            DisplayName=['$\varepsilon = 24 \rm{m^3/s^2}, d = ',num2str(0.15*i),'\rm{mm}$'] );
        hold on
        % 记录图例标签（示例格式）
        % legendEntries{idx} = sprintf('A%d, B%d', iA, iB);

end



% alpha = 2;      % 替换为所需的alpha值
% beta = 20;     % 替换为所需的beta值
% ratio = 0;
% f_v = linspace(0, 1, 1000);  % 关键修改
% pdf = (1-ratio) * (betapdf(f_v, alpha, beta) + betapdf(1 - f_v, alpha, beta));
% plot(f_v, pdf, 'g:', 'LineWidth', 2, 'DisplayName', 'Approximate PDF');
legend(Interpreter='latex');
xlabel('$f_v$',Interpreter='latex')
ylabel('$\beta(f_v;d,\varepsilon)$',Interpreter='latex')
xlim([0 0.5]);
ylim([0 10])
daspect([1/10 1/0.5 1])
grid on
% set(gca,'YScale','log')



%%
ratio = 0;
alpha = 2;      % 替换为所需的alpha值
beta = 40;     % 替换为所需的beta值
sigma = 0;
f_v = linspace(0, 1, 1000);  % 关键修改
pdf = (1-ratio) * (betapdf(f_v, alpha, beta) + betapdf(1 - f_v, alpha, beta))/2 ...
      + ratio * unifpdf(f_v, sigma, 1-sigma);

plot(f_v, pdf, 'k--', 'LineWidth', 1, 'DisplayName', 'Bimodal DDSD');
hold on

ratio = 0.2;
pdf = (1-ratio) * (betapdf(f_v, alpha, beta) + betapdf(1 - f_v, alpha, beta))/2 ...
      + ratio * normpdf(f_v, 0.5, 0.3);

plot(f_v, pdf, 'k-.', 'LineWidth', 1, 'DisplayName', 'Bimodal DDSD with non-zero probability of equal breakup');


plot_ao_pdf(0.15, 0.2);


hold off;
grid on;

% % 添加图例（可选）
legend('Location', 'northeast', 'Interpreter','latex');
text(0,4.8,"(a)",Interpreter="latex",FontSize=20)
xlabel('$f_v = (d/d_0)^3$',Interpreter='latex')
ylabel('$\beta(f_v;d,\varepsilon)$',Interpreter='latex')
legend(Interpreter='latex')
ylim([0,5])
daspect([1 5 1])

function plot_ao_pdf(a, h_low)
% 绘制凹字形分布的理论PDF曲线
% 输入:
%   a     : 左侧区间结束点 (0 < a < 0.5)
%   h_low : 中间区间(a,1-a)的概率密度

% 参数验证
if a <= 0 || a >= 0.5
    error('a必须在(0,0.5)区间内');
end

% 计算对称点b和两侧密度h_high
b = 1 - a;
h_high = (1 - (1 - 2*a)*h_low) / (2*a);

% 检查密度有效性
if h_high <= 0
    error('h_low值过大导致h_high<=0! 最大允许h_low=%.4f', 1/(1-2*a));
end

% 创建PDF的关键点 (避免跳跃处的垂直线)
x = [0, a, a, b, b, 1];  % x坐标
y = [h_high, h_high, h_low, h_low, h_high, h_high];  % y坐标

% 绘制PDF曲线
% figure;
plot(x, y, 'k:', 'LineWidth', 1, 'DisplayName', 'Mixed Uniform DDSD');
grid on;

% 添加标记点和参考线
% hold on;
% plot([a, a], [0, max(h_high, h_low)*1.1], 'r--', 'LineWidth', 0.5); % a处参考线
% plot([b, b], [0, max(h_high, h_low)*1.1], 'r--', 'LineWidth', 0.5); % b处参考线
% plot(a, h_high, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % 左跳跃点
% plot(a, h_low, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');  % 左跳跃点
% plot(b, h_low, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');  % 右跳跃点
% plot(b, h_high, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % 右跳跃点
% 
% 设置坐标轴和标签
% axis([0 1 0 max(h_high, h_low)*1.1]);
% xlabel('x');
% ylabel('Probability Density');
% title(sprintf('理论PDF曲线 (a=%.2f, h_{low}=%.2f)', a, h_low));
% legend('PDF', '分段点', 'Location', 'best');

% 添加密度值标注
% text(0.1*a, h_high, sprintf('h_{high}=%.2f', h_high), 'FontSize', 10);
% text((a+b)/2, h_low, sprintf('h_{low}=%.2f', h_low), 'FontSize', 10);
% hold off;
end

function smoothedData = smoothCurve(data)
    % 输入：data - N×2矩阵，第一列为横坐标x，第二列为纵坐标y
    % 输出：smoothedData - 平滑后的数据，格式同输入
    
    % 确保数据按横坐标排序
    [~, idx] = sort(data(:,1));
    sortedData = data(idx, :);
    x = sortedData(:,1);
    y = sortedData(:,2);
    
    % 自动计算最优窗口长度（奇数）
    n = length(x);
    windowSize = min(floor(n/3)*2 + 1, 15); % 不超过15的奇数
    
    % 使用Savitzky-Golay滤波器（多项式阶数=2）
    if windowSize > 2
        smoothedY = sgolayfilt(y, 2, windowSize);
    else
        warning('数据点过少，无法有效平滑');
        smoothedY = y;
    end
    
    % 返回平滑后数据
    smoothedData = [x, smoothedY];
end