   function tau = calculate_tau(t_c, d)
    % 给定参数
    mu_d = 0.0021;%2.4e-6*860;     
    sigma = 0.004;    
    theta_c = 1;  

    if(0) % 简单模型，没有非线性项
        t_ca = d * mu_d /sigma;% 张力特征时间
    
        % 初始化 tau 数组
        tau = zeros(size(t_c));
    
        % 计算 tau 的值
        for i = 1:length(t_c)
            % 根据公式计算 tau
            exponent = exp( - t_c(i)/t_ca);
            tau(i) = theta_c * mu_d/(t_ca*(1 - exponent));
        end
    else

        syms x(t) a c k
        eqn = k*diff(x,t) + a*x*(1 - x) == c;
        sol = dsolve(eqn,x(0)==0);
        res = solve(sol==1);
        rs2 = subs(res,a,sigma/d);
        rs3 = subs(rs2,k,mu_d);
        parfor i = 1:length(t_c)
            tau(i) = vpasolve(rs3==t_c(i));
                    % disp(i)
        end

    end
end