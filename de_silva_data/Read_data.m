close all
clear all
clc
load('AEM_data_VEL_FIELDS_progress_save2_Rahul_0.12.mat');

U_ALL = load_data.U_ALL;
V_ALL = load_data.V_ALL/utau_est;
W_ALL = load_data.W_ALL/utau_est;
Num_layers = 1:6;

mean_U = squeeze(mean(mean(sum(U_ALL,4)/utau_est+Uinf,2),1));
zpos = (0:size(U_ALL,3)-1)*32; % grid spacing is 32 wall units
zpos_delta = zpos/3200;
 
xpos = (0:size(U_ALL,2)-1)*32;
xpos_delta = xpos./3200;

ypos = (0:size(U_ALL,1)-1)*32;
ypos_delta = ypos./3200;

%% add up all hier. to get full field
close all
U = sum(U_ALL,4)/utau_est+Uinf;
V = sum(V_ALL,4)/utau_est;
W = sum(W_ALL,4)/utau_est;
save("example_data.mat","U","V","W","zpos_delta","xpos_delta","ypos_delta")%,"utau_est","Uinf")

% Test plots
figure;imagesc(xpos_delta,zpos_delta,squeeze(U(1,:,:))'); daspect([1 1 1]);set(gca,'ydir','normal') % Streamwise - wall normal plane (rotated)
figure;imagesc(xpos_delta,ypos_delta,squeeze(U(:,:,5))'); daspect([1 1 1]); %  wall parallel plane 

%%
figure;

for i = 1:101
    mean_velocity(i) = mean(mean(U(:,:,i)));
end
% mean_velocity = mean_velocity - mean_velocity(1)*ones(1,101);

semilogx(zpos,mean_velocity,'-x',linewidth=2,DisplayName='Raw Data');
hold on 
semilogx(zpos,log(zpos)/0.41+7.5,DisplayName='Log Law');
xlabel('Wall-normal Distance $z^+$',Interpreter='latex');
ylabel('Streamwise Velocity $U^+$',Interpreter='latex');
legend();
hold off