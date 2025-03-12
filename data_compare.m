%检验数据格式是否一致
load("example_data.mat")

% Test plots
figure;imagesc(xpos_delta,zpos_delta,squeeze(U(1,:,:))'); daspect([1 1 1]);set(gca,'ydir','normal') % Streamwise - wall normal plane (rotated)
figure;imagesc(xpos_delta,ypos_delta,squeeze(U(:,:,5))'); daspect([1 1 1]); %  wall parallel plane 
%x是流向，y是展向，z是垂向

%%
clear,clc;
load("self_example_data.mat")

% Test plots
figure;imagesc(xpos_delta,zpos_delta,squeeze(U(1,:,:))'); daspect([1 1 1]);set(gca,'ydir','normal');colorbar; % Streamwise - wall normal plane (rotated)
figure;imagesc(xpos_delta,ypos_delta,squeeze(U(:,:,5))); daspect([1 1 1]);colorbar; %  wall parallel plane 

