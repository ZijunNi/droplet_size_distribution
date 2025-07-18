% MATLAB script to compute daughter bubble/droplet size distribution as per Wang et al. (2003)
% Chemical Engineering Science 58 (2003) 4629-4637
% Fixed normalization and small droplet size handling

% Parameters
d = 0.001; % Mother bubble/droplet diameter (m), e.g., 0.1 mm
sigma = 0.00725; % Surface tension (N/m) for air-water system
epsilon = 6.0; % Energy dissipation rate (m^2/s^3)
rho_c = 860; % Density of continuous phase (kg/m^3)
alpha_d = 0.01; % Volume fraction of dispersed phase
n = 1; % Bubble/droplet number density (1/m^3), assumed for simplicity
delta = 0.001; % Model parameter to avoid singularity
nu = 2.4e-6; % Kinematic viscosity of water (m^2/s)
lambda_d = (nu^3 / epsilon)^0.25; % Kolmogorov length scale
lambda_min_ratio = 11.4; % lambda_min / lambda_d
lambda_min = lambda_d * lambda_min_ratio; % Minimum eddy size (m)

% Check maximum stable diameter
d_max = (sigma / (rho_c * epsilon^(2/3)))^(3/5) * 0.725; % Approximate d_max, constant from Hinze (1955)
if d < d_max
    warning('Droplet diameter (d = %.2e m) is smaller than maximum stable diameter (d_max = %.2e m). Breakup probability is very low.', d, d_max);
end

% Breakup fraction range
fv = linspace(0.001, 0.5, 100); % Breakup fraction, avoid f_v = 0 to prevent singularity
b_fv = zeros(size(fv)); % Initialize breakup rate b(f_v|d)

% Check if d < lambda_min
if d < lambda_min
    warning('Droplet diameter (d = %.2e m) is smaller than minimum eddy size (lambda_min = %.2e m). No breakup occurs.', d, lambda_min);
    beta_fv = zeros(size(fv)); % Return zero distribution
else
    % Calculate breakup rate b(f_v|d) for each f_v
    for i = 1:length(fv)
        f_v = fv(i);
        % Numerical integration over lambda from lambda_min to d
        lambda = linspace(lambda_min, d, 1000);
        integrand = zeros(size(lambda));
        
        for j = 1:length(lambda)
            lambda_j = lambda(j);
            
            % Turbulent velocity of eddy (Eq. 3)
            u_lambda = sqrt(2) * (epsilon * lambda_j)^(1/3);
            
            % Eddy number density (Eq. 4)
            n_lambda = 0.822 * (1 - alpha_d) / lambda_j^4;
            
            % Collision frequency density (Eq. 5)
            omega_lambda = 0.923 * (1 - alpha_d) * n * epsilon^(1/3) * (lambda_j + d)^2 / lambda_j^(11/3);
            
            % Mean eddy kinetic energy (Eq. 13)
            e_bar_lambda = (pi/6) * lambda_j^3 * rho_c * (u_lambda^2 / 2);
            
            % Energy distribution integration for P_b(f_v|d,lambda) (Eq. 14)
            e_lambda = linspace(0, 100*e_bar_lambda, 200); % Extended energy range for small droplets
            P_b = 0;
            for k = 1:length(e_lambda)
                e = e_lambda(k);
                
                % Minimum breakup fraction (Eq. 7)
                f_v_min = ((pi * lambda_j^3 * sigma) / (6 * e * d))^3;
                if f_v_min > 0.5
                    f_v_min = 0; % Ensure physically meaningful range
                end
                
                % Surface energy increase (Eq. 8)
                c_f = f_v^(2/3) + (1 - f_v)^(2/3) - 1;
                Delta_e_i = c_f * pi * d^2 * sigma;
                
                % Maximum breakup fraction (Eq. 10)
                c_f_max = min((2^(1/3) - 1), e / (pi * d^2 * sigma));
                if c_f_max <= 0
                    f_v_max = 0;
                else
                    % Solve for f_v_max where c_f(f_v_max) = c_f_max
                    try
                        f_v_max = fzero(@(x) x^(2/3) + (1-x)^(2/3) - 1 - c_f_max, [0.001, 0.5]);
                    catch
                        f_v_max = 0; % Handle numerical issues
                    end
                end
                
                % Breakup probability (Eq. 11)
                if f_v_max - f_v_min >= delta && f_v >= f_v_min && f_v <= f_v_max
                    P_b_fv_e_lambda = 1 / (f_v_max - f_v_min);
                else
                    P_b_fv_e_lambda = 0;
                end
                
                % Energy distribution (Eq. 12)
                P_e = (1 / e_bar_lambda) * exp(-e / e_bar_lambda);
                
                % Accumulate P_b(f_v|d,lambda)
                P_b = P_b + P_b_fv_e_lambda * P_e * (e_lambda(2) - e_lambda(1));
            end
            
            % Integrand for b(f_v|d)
            integrand(j) = P_b * omega_lambda;
        end
        
        % Numerical integration over lambda (Eq. 15)
        b_fv(i) = trapz(lambda, integrand);
    end
    
    % Total breakup rate (Eq. 16, adjusted for symmetry)
    b_total = trapz(fv, b_fv);
    
    % Daughter size distribution (Eq. 17, fixed normalization)
    if b_total > 0
        beta_fv = b_fv / b_total; % Normalized over [0, 0.5]
    else
        beta_fv = zeros(size(fv)); % Handle case with no breakup
    end
end

% Plot the daughter size distribution
figure;
plot(fv, beta_fv, 'b-', 'LineWidth', 2);
xlabel('Breakup Fraction f_v');
ylabel('Daughter Size Distribution \beta(f_v, d)');
title(sprintf('Daughter Size Distribution (d = %.2e m, \\epsilon = %.1f m^2/s^3)', d, epsilon));
grid on;

% Save the results to a file
results = [fv; beta_fv]';
save('daughter_size_distribution_fixed.txt', 'results', '-ascii');