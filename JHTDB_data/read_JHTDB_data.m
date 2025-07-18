clear,clc

num_datasets = 1;%需要读取的数据集数量

for i = 1:num_datasets
    file_name = ['JHTDB_Re_tau_1000_test_',num2str(i),'.mat'];
res{i} = load(file_name);
U_read{i} = res{i}.results;
end

z_point = [];
U = [];
for i = 1:num_datasets
    z_point = [z_point res{i}.z_points];
    U = cat(1, U, U_read{i});
end
x_point = res{1}.x_points;
y_point = res{1}.y_points;

figure;imagesc(x_point,z_point,squeeze(U(:,10,:))'); %daspect([1 1 1]); %  wall parallel plane 

%%

U = permute(U, [1 3 2]);
V = zeros(size(U));
W = zeros(size(U));
xpos_delta = x_point;
ypos_delta = z_point;
zpos_delta = y_point;
%%
% Test plots
% figure;imagesc(xpos_delta,zpos_delta,squeeze(U(1,:,:))'); %daspect([1 1 1]);set(gca,'ydir','normal') % Streamwise - wall normal plane (rotated)
figure;imagesc(xpos_delta,ypos_delta,squeeze(U(:,:,1))'); %daspect([1 1 1]); %  wall parallel plane 
    mean_U = squeeze(mean(mean(U,2),1));%save mean velocity profile
    plot(zpos_delta*1000,mean_U)


%%

save("data_Re_tau_1000_JHTDB.mat","U","V","W","xpos_delta","ypos_delta","zpos_delta");