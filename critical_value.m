function [critical_value,physical_duration,physical_threshold,physical_tau,alternative_density] = ...
    critical_value(U,V,dx,zpos_delta,Reynolds_number,ratio,left_bound,right_bound)

    rho_c = [1,1000];%mass density of comtinumm phase

    %%%% density = 1 %%%%
        myFunction = @(x) condition_function(U,V,dx,zpos_delta,Reynolds_number,ratio,x,rho_c(1)); 
        tolerance = 1e-9; 
        left_bound0 = left_bound;
        right_bound0 = right_bound;

        if myFunction(left_bound) == myFunction(right_bound)
            error('Initial Bound Error 1, [left_bound, right_bound] need to be adjusted');
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
        fprintf('The Critical Value is %.8f\n', critical_value);


    %%%% density = 1000 %%%%
        myFunction = @(x) condition_function(U,V,dx,zpos_delta,Reynolds_number,ratio,x,rho_c(2)); 
        tolerance = 1e-9; 
        left_bound = left_bound0;
        right_bound = right_bound0;
        
        if myFunction(left_bound) == myFunction(right_bound)
        error('Initial Bound Error 2, [left_bound, right_bound] need to be adjusted');
        end
        
        while (right_bound - left_bound) > tolerance
            mid = (left_bound + right_bound) / 2; % 取中点
            if myFunction(mid) == myFunction(left_bound)
                left_bound = mid; % 如果中点值与左端点值相同，更新左端点
            else
                right_bound = mid; % 否则更新右端点
            end
        end
        
        alternative_density{1} = (left_bound + right_bound) / 2;
        [~,alternative_density{2},alternative_density{3},alternative_density{4}] = condition_function(U,V,dx,zpos_delta,Reynolds_number,ratio,critical_value,rho_c(2));


end