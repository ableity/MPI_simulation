function [out_x,out_y]=MPI_sim_2D_2(img)
%采样频率
f_s_raw = 1250000;
%采样时间
t = 1/2500;
%激励频率
f = 25000;
% 笛卡尔轨迹N（倍率）
N = 10;

%外加磁场峰值，单位mT
H_peak_x = 12.5;
H_peak_y = 12.5;
%外加梯度磁场梯度,单位T/m
H_gradient = 2.5;
%峰值单位转化
H_peak_x = H_peak_x*1e-3;
H_peak_y = H_peak_y*1e-3;
%外加磁场频率
fx = f;
fy = f/N;
% 时间序列
% 采样率放大了10倍，采样时间放大了3倍
f_s = 10*f_s_raw;
t_real = 1/f_s:1/f_s:t*3;
% time_pixel为时间细分数，本程序中时间的最小单位
[~,time_pixel] = size(t_real);

%FOV细分数量，即单个方向像素个数
num_of_pixel = 100; 
%计算FOV，绝对零磁场点最大移动范围再扩大1.2倍（或许这里有点小）
FOV = max(H_peak_x,H_peak_y)/H_gradient*2;
FOV = FOV*1.2;
%像素尺寸
pixel = FOV/num_of_pixel;
%梯度场空间分布，
G_x= -0.5*H_gradient*FOV+pixel*H_gradient:pixel*H_gradient:0.5*H_gradient*FOV;
G_y= -0.5*H_gradient*FOV+pixel*H_gradient:pixel*H_gradient:0.5*H_gradient*FOV;
%激励场时间分布
D_x = H_peak_x*sin(2*pi*fx*t_real);
D_y = H_peak_y*sin(2*pi*fy*t_real);
% 建立两个矩阵，一行数据是在相同时间不同位置下的H
% 一列数据是在相同位置下不同时间的t
G_x =  repmat(G_x,1,100);
G_x = reshape(G_x,1,[]);
G_x =  repmat(G_x,time_pixel,1);

G_y =  repmat(G_y,100,1);
G_y = reshape(G_y,1,[]);
G_y =  repmat(G_y,time_pixel,1);

D_x = reshape(D_x,[],1);
D_x =  repmat(D_x,1,num_of_pixel^2);

D_y = reshape(D_y,[],1);
D_y =  repmat(D_y,1,num_of_pixel^2);

% 梯度场、激励场合成
H_x = D_x+G_x;
H_y = D_y+G_y;

%仿体，img外部输入
phantom = reshape(img,1,[]);
phantom_real = repmat(phantom,time_pixel,1);
%仿真
[out_x,out_y] = MPI_sim_2D_1(H_x,H_y,f_s,phantom_real);


% 截取一个周期并降采样，前文将采样时间放大了3倍，这里取最中间的周期作为输出
out_x = out_x(round(t*f_s+1):round(2*t*f_s));
% 下采样10倍，和前文的上采样10倍对应
out_x= downsample(out_x,10);

% 同x通道
out_y = out_y(round(t*f_s+1):round(2*t*f_s));
out_y= downsample(out_y,10);


% x = 1:100;
% x = repmat(x,1,100);
% y = 1:100;
% y = repmat(y,100,1);
% y = reshape(y,1,10000);
% u = repmat(G_x,1,100);
% v = repmat(G_y,100,1);
% v = reshape(v,1,10000);

% figure
% for i = 1:15000
%     quiver(x,y,u+D_x(i),v+D_y(i))
%     pause(0.001)
% end

