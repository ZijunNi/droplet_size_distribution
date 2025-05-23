function [critical_value,physical_duration,physical_threshold,physical_tau] = ...
    critical_value(U,V,dx,zpos_delta,Reynolds_number,ratio,left_bound,right_bound)

        rho_c = 1000;%mass density of comtinumm phase

        myFunction = @(x) condition_function(U,V,dx,zpos_delta,Reynolds_number,ratio,x,rho_c); 
        tolerance = 1e-9; 

        if myFunction(left_bound) == myFunction(right_bound)
            error('Initial Bound Error, [left_bound, right_bound] need to be adjusted');
        end
        
        while (right_bound - left_bound) > tolerance
            mid = (left_bound + right_bound) / 2; % 取中点
            if myFunction(mid) == myFunction(left_bound)
                left_bound = mid; % 如果中点值与左端点值相同，更新左端点
            else
                right_bound = mid; % 否则更新右端点
            end
        end
        
        critical_value = (left_bound + right_bound) / 2;
        [~,physical_duration,physical_threshold,physical_tau] = condition_function(U,V,dx,zpos_delta,Reynolds_number,ratio,critical_value,rho_c(1));
        fprintf('The Critical Value is %6.5e\n', critical_value);

end