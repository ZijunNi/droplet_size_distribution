% MATLAB code to calculate daughter particle size distribution eta(v:vf) for air-water system
% Based on the model from Luo and Svendsen (1996), corrected for normalization

clear all;
close all;

% Physical parameters for air-water system (SI units)
d = 0.0001; % Parent particle diameter (m), e.g., 3 mm
epsilon = 5; % Energy dissipation rate (m^2/s^3)
sigma = 0.0052; % Surface tension (N/m) for air-water
rho_c = 860; % Continuous phase density (kg/m^3)
beta = 2.0; % Turbulence constant (Kuboi et al., 1972a, more accurate value)
c4 = 0.923; % Constant c4 (Eq. 18, c4 = c3 * pi * beta^0.5 / 4)
epsilon_d = 0; % Dispersed phase volume fraction (assumed negligible)
xi_min = 0.2; % Minimum eddy size ratio (lambda_min/d)
n = 1; % Number density of particles (m^-3, assumed for simplicity)

% Breakage volume fraction f_BV range
f_BV = linspace(0.01, 0.99, 100); % Avoid singularities at 0 and 1
eta_v = zeros(size(f_BV)); % Initialize eta*v distribution

% Calculate partial breakage rate Omega_B(v:vf_BV)
for i = 1:length(f_BV)
    % Surface area increase coefficient c_f
    cf = f_BV(i)^(2/3) + (1 - f_BV(i))^(2/3) - 1;
    
    % Integrand for Omega_B
    integrand = @(xi) (1 + xi).^2 .* xi.^(-11/3) .* ...
        exp(-12 * cf * sigma ./ (beta * rho_c * epsilon^(2/3) * d^(5/3) * xi.^(11/3)));
    Omega_B = c4 * (1 - epsilon_d) * n * (epsilon / d^2)^(1/3) * ...
        integral(integrand, xi_min, 1, 'AbsTol', 1e-8, 'RelTol', 1e-8);
    
    % Store for numerator of eta
    eta_v(i) = Omega_B;
end

% Calculate total breakage rate Omega_B(v) for normalization
total_integrand = @(xi, f) (1 + xi).^2 .* xi.^(-11/3) .* ...
    exp(-12 * (f.^(2/3) + ( one_minus_f(f) ).^(2/3) - 1) * sigma ./ ...
    (beta * rho_c * epsilon^(2/3) * d^(5/3) * xi.^(11/3)));
denominator = integral2(total_integrand, xi_min, 1, 0, 0.5, 'AbsTol', 1e-8, 'RelTol', 1e-8);
Omega_B_total = c4 * (1 - epsilon_d) * n * (epsilon / d^2)^(1/3) * denominator;

% Calculate eta(v:vf) * v
v = (pi / 6) * d^3; % Volume of parent particle
eta_v = eta_v / Omega_B_total; % Normalize by total breakage rate

% Plot the daughter particle size distribution
% figure;
plot(f_BV, eta_v, 'b-', 'LineWidth', 2);
xlabel('Breakage Volume Fraction f_{BV}');
ylabel('\eta(v:v_f) \cdot v');
title('Daughter Particle Size Distribution (Air-Water, d=3mm, \epsilon=0.5 m^2/s^3)');
grid on;
hold on

% Optional: Compare with d=6mm, epsilon=1 m^2/s^3
% d2 = 0.006;
% epsilon2 = 1;
% eta_v2 = zeros(size(f_BV));
% for i = 1:length(f_BV)
%     cf = f_BV(i)^(2/3) + (1 - f_BV(i))^(2/3) - 1;
%     integrand = @(xi) (1 + xi).^2 .* xi.^(-11/3) .* ...
%         exp(-12 * cf * sigma ./ (beta * rho_c * epsilon2^(2/3) * d2^(5/3) * xi.^(11/3)));
%     Omega_B = c4 * (1 - epsilon_d) * n * (epsilon2 / d2^2)^(1/3) * ...
%         integral(integrand, xi_min, 1, 'AbsTol', 1e-8, 'RelTol', 1e-8);
%     eta_v2(i) = Omega_B;
% end
% denominator2 = integral2(total_integrand, xi_min, 1, 0, 0.5, 'AbsTol', 1e-8, 'RelTol', 1e-8);
% Omega_B_total2 = c4 * (1 - epsilon_d) * n * (epsilon2 / d2^2)^(1/3) * denominator2;
% v2 = (pi / 6) * d2^3;
% eta_v2 = eta_v2 / Omega_B_total2;
% 
% % Add second curve to plot
% hold on;
% plot(f_BV, eta_v2, 'r--', 'LineWidth', 2);
% legend('d=3mm, \epsilon=0.5 m^2/s^3', 'd=6mm, \epsilon=1 m^2/s^3');
% hold off;

% Helper function for 1-f in integrand
function y = one_minus_f(f)
    y = 1 - f;
end
