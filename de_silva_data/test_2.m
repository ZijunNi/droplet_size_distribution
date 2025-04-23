z_plus = squeeze(z_frame(1,1,:)) * Re_tau;
logical_region = (z_plus > 50) & (z_plus < max(0.1 * Re_tau,60));
indices = find(logical_region);
p_lin = polyfit(log(z_plus(indices)), squeeze(u_mean(indices)), 1);
% p_lin = polyfit(log(squeeze(z_frame(1,1,begin_pos:begin_pos+3))*Re_tau),squeeze(u_mean(begin_pos:begin_pos+3)),1);  % 2:5 corresponds to the range in the mean flow profile I fit too
figure;semilogx(squeeze(z_frame(1,1,:))*Re_tau,squeeze(u_mean),'-x');hold on;
plot([1 Re_tau],polyval(p_lin,log([1 Re_tau])),'--k');
plot([1 Re_tau],polyval([2.5 p_lin(2)],log([1 Re_tau])),'--b');
Uinf_est = 5-p_lin(2);
% get Uinf estimates
Uinf = polyval([2.5 5],log(Re_tau));
U_ALL = U_ALL +Uinf;

%%
% add up all hier. to get full field
close all

U = sum(U_ALL(:,:,:,1:n),4)/utau_est;%+Uinf;
V = sum(V_ALL(:,:,:,1:n),4)/utau_est;
W = sum(W_ALL(:,:,:,1:n),4)/utau_est;