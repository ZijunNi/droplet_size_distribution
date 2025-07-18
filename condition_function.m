function [result,physical_duration,physical_threshold,physical_tau] = ...
condition_function(u,v,dx,zpos_delta,Reynolds_number,ratio,diameter,rho_c)
    
        delta = 0.5e-2;%半槽宽
        nu_c = 1.8e-6;% 水的运动粘度nu
        u_tau = Reynolds_number*nu_c/delta;
        
        [threshold_series, duration] = threshold_line(u, v, dx, ratio, diameter, zpos_delta, Reynolds_number);

        physical_threshold = threshold_series*0.5*rho_c*u_tau^2;%计算有量纲的阈值 
        
        %计算得到该粒子半径对应的阈值-时间线
        %-----使用当地平均速度作为u_bar
                % interped_2D_field = u_tau*extract_2d_slice_x_interp(u,zpos_delta,diameter,Reynolds_number);%导出考察位置处的速度切片
                % u_bar = mean(interped_2D_field(:));
        %-----使用0.8*U_inf作为u_bar
                mid_line_loc = find(zpos_delta-1<=0,1,'last');%中线位置
                mid_plane = squeeze(u(:,:,mid_line_loc));
                if(Reynolds_number == 180)
                    u_bar = mean(mid_plane(:));
                elseif(Reynolds_number == 200)
                    u_bar = 18.534401289144565;% From CM de Silva dataset
                elseif(Reynolds_number == 400)
                    u_bar = 19.978661;% From CM de Silva dataset
                elseif(Reynolds_number == 800)
                    u_bar = 21.711529;% From CM de Silva dataset
                elseif(Reynolds_number == 1600)
                    u_bar = 23.444397;% From CM de Silva dataset
                elseif(Reynolds_number == 1000)
                    u_bar = 22.638488632724943;% JHTDB dataset
                end
                u_bar = 0.5*u_bar;
        %-----使用9.5*u_tau作为u_bar
                % u_bar =9.50*u_tau;

        physical_duration = duration.*delta/u_bar;%用当地流向速度将空间坐标转化为时间坐标
        physical_tau = calculate_tau(physical_duration,delta*diameter);% 计算破碎线
        
        % figure;
        % disp(u_bar)
        % plot(physical_duration,physical_tau,"DisplayName",'Breakup line')
        % hold on
        % plot(physical_duration,physical_threshold,"DisplayName",'Threshold line')
        % hold off

        result = has_element_less_than(physical_tau,physical_threshold);
end
