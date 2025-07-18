% 本程序处理Re_tau向Re的转化
clear,clc

K = 0.74;
r_i = 25e-3;
d = 10e-3;

Re_tau = [200 400 800 1000 1600];

Re = ((r_i/(0.5*d*sqrt(K)))*Re_tau).^(1/0.79);
