%将自有程序运行数据转化为液滴碎裂程序可用的格式
clc,clear;
load("fluid_data.mat")
%x应为流向，y是展向，z是垂向

xpos_delta = x;
ypos_delta = z;
zpos_delta = y_grid;
Re_tau = 180;
Re_bulk = 2820;
u_tau = Re_tau/Re_bulk;


U = zeros(128,128,128*length(data_xz_fluid_U));
V = zeros(128,128,128*length(data_xz_fluid_U));
W = zeros(128,128,128*length(data_xz_fluid_U));
for i = 1:length(data_xz_fluid_U)
    U(:,:,128*(i-1)+1:128*i) = data_xz_fluid_U{i};
    V(:,:,128*(i-1)+1:128*i) = data_xz_fluid_V{i};
    W(:,:,128*(i-1)+1:128*i) = data_xz_fluid_W{i};
end

U = permute(U,[3,1,2])/u_tau;
V = permute(V,[3,1,2])/u_tau;
W = permute(W,[3,1,2])/u_tau;

mid = V;
V = W;
W = mid;

save("self_example_data.mat","U","V","W","zpos_delta","xpos_delta","ypos_delta")


