function sci_str = num2sci(num, n_significant)
% 将小数转换为科学计数法字符串（TeX格式）
% 输入：
%   num: 要转换的数字（标量）
%   n_significant: 有效数字位数（默认3位）
% 输出：
%   sci_str: 科学计数法字符串（如 '1.23 \times 10^{2}'）

    if nargin < 2
        n_significant = 3; % 默认有效数字位数
    end

    if num == 0
        sci_str = '0';
        return;
    end

    % 计算指数和底数
    exponent = floor(log10(abs(num)));
    mantissa = num / 10^exponent;

    % 四舍五入到指定有效数字
    factor = 10^(n_significant - 1);
    mantissa = round(mantissa * factor) / factor;

    % 调整底数范围（确保 |mantissa| ∈ [1, 10)）
    if abs(mantissa) >= 10
        mantissa = mantissa / 10;
        exponent = exponent + 1;
    end

    % 格式化底数字符串
    if n_significant == 1
        fmt = '%.0f'; % 无小数位
    else
        fmt = ['%.' num2str(n_significant - 1) 'f']; % 保留n-1位小数
    end
    mantissa_str = sprintf(fmt, mantissa);

    % 组合成科学计数法字符串（TeX格式）
    if exponent == 0
        sci_str = mantissa_str;
    else
        sci_str = [mantissa_str ' \times 10^{' num2str(exponent) '}'];
    end
end