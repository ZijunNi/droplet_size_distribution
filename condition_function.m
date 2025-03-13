function [result,physical_duration,physical_threshold,physical_tau] = ...
condition_function(u,v,dx,zpos_delta,Reynolds_number,ratio,diameter,rho_c)
    
        delta = 0.5e-2;%半槽宽
        nu_c = 1.8e-6;% 水的运动粘度nu
        u_tau = Reynolds_number*nu_c/delta;
        
        [threshold_series, duration] = threshold_line(u, v, dx, ratio, diameter,zpos_delta,Reynolds_number);
        physical_threshold = threshold_series*0.5*rho_c*u_tau^2;%计算有量纲的阈值
        
        %计算得到该粒子半径对应的阈值-时间线
        interped_2D_field = u_tau*extract_2d_slice_x_interp(u,zpos_delta,diameter,Reynolds_number);%导出考察位置处的速度切片
        physical_duration = duration.*delta/mean(interped_2D_field(:));%用当地流向速度将空间坐标转化为时间坐标
        physical_tau = calculate_tau(physical_duration,delta*diameter);% 计算破碎线
        
        % semilogx(physical_duration,physical_threshold,'-s');
        % hold on 
        % semilogx(physical_duration,physical_tau,'-x');
        % t = ['$d/\delta=$ ',num2str(diameter,'%6.5e')];
        % title(t,'interpreter','latex')
        % hold off
        
        result = has_element_less_than(physical_tau,physical_threshold);
end
