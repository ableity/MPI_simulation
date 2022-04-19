clc
clear
close all

%仿体
phantom = zeros(1,1000);
phantom(50)= 1;
%仿真
out = MPI_sim_1D_2(phantom);
figure
plot(out)
title("接收的信号")