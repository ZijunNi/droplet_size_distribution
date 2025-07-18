% 自定义聚合概率函数 
function [p,K] = aggregation_prob(s1, s2, break_threshold, Re_tau)

    %%%%%%% 实验参数 %%%%%%%
        delta = 0.5e-2;
        nu_c = 1.8e-6;
        rho_c = 1000;
        rho_d = 866;
        sigma = 0.004;
    %%%%%%% 实验参数 %%%%%%%

    %%%%%%% 模型参数 %%%%%%%
        K_1 = 1.2e-4;%2e-4
        K_2 = 0;%0.5e-6;
        K_3 = 0;%5e-3;
        K = [K_1,K_2,K_3];
    %%%%%%% 模型参数 %%%%%%%

    v1 = (s1/break_threshold)^3;
    v2 = (s2/break_threshold)^3;
    % 基础概率与尺寸乘积成正比，保证小液滴聚合概率低
    % p = min(0.0005, 0.01*(s1/breakiyou_threshold)*(s2/break_threshold));2358.mat
    dissipation = Re_tau^4*nu_c^3/delta^4;


    collision = K_1*(v1^(2/3) + v2^(2/3))*sqrt(v1^(2/9)+v2^(2/9))*dissipation^(1/3);
    efficiency = exp(-K_2*nu_c*rho_c^2*dissipation*(v1^(1/3)*v2^(1/3)/(v1^(1/3)+v2^(1/3)))^4/(sigma^2));
    
    v = sqrt(v1*v2);
    breakup_rate = v^(-2/9)*dissipation^(1/3)*exp(-K_3*sigma/(rho_d*dissipation^(2/3)*v^(5/9)));

    if(K_2==0&&K_3==0)
        p = collision;
    elseif(K_3==0)
        p = collision*efficiency;
    else
        p = collision*efficiency/breakup_rate;
    end
    if(p>0.99)
        warning('最大概率接近1，无法形成液滴群。')
    end
    %Coulaloglou C A, Tavlarides L L. Description of interaction processes in agitated liquid-liquid dispersions[J]. Chemical Engineering Science, 1977, 32(11): 1289-1297.

end