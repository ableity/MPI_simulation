clc
clear
close all
%二维MPI仿真脚本
%定义图像仿体（尺寸只能是100*100）
img = zeros(100,100);
img(20,20) = 10;
img(50,50) = 10;
%仿真
tic
[out_x,out_y] = MPI_sim_2D_2(img);
toc

%显示信号
figure
subplot(3,4,[1 2 5 6])
imshow(img,[])
title("仿体")
subplot(3,4,[9 10])
plot(out_x)
title("X通道信号")
subplot(3,4,[11 12])
plot(out_y)
title("y通道信号")


t_real = 1/12500000:1/12500000:4.e-04;
D_x = 100*sin(2*pi*25000*t_real);
D_y = 100*sin(2*pi*2500*t_real);
D_x= downsample(D_x,10);
D_y= downsample(D_y,10);
subplot(3,4,[3 4 7 8])
comet(D_x,D_y)
title("FFP轨迹")






