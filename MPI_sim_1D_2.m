function out = MPI_sim_1D_2(phantom)
% MPI一维仿真（二次封装）
%% 根据外磁场和仿体总结参数
% 本部分为外场所有参数，此部分可以单独写脚本当外置输入，所以分开写
H_peak = 12.5;
%外加梯度磁场梯度,单位T/m
H_gradient = 2.5;
%峰值单位转化
H_peak = H_peak*1e-3;
%采样频率
f_s_raw = 1250000;
%采样时间
t = 1/25000;
%激励频率
f = 25000;
%FOV细分数量(决定了体素大小)
num_of_pixel = 1000; 
%FOV,一个长度，单位为m
FOV = H_peak/H_gradient*2;
FOV = FOV*1.2;
%像素长度 
pixel = FOV/num_of_pixel;
% 增加采样率与采样点数
% 时间序列（按照采样率的10倍，后面通过降采样回到真实采样率）
% 时间用了3倍，后面可以取中间的一段

% 激励磁场
f_s = 10*f_s_raw;
t_real = 1/f_s:1/f_s:t*3;
% 外加激励磁场随时间的变化
H_t = H_peak*sin(2*pi*f*t_real);
% 梯度磁场
H_x = -0.5*H_gradient*FOV+pixel*H_gradient:pixel*H_gradient:0.5*H_gradient*FOV;


% 建立一个矩阵，一行数据是在相同时间不同位置下的H
% 一列数据是在相同位置下不同时间的t
[~,H_t_SIZE] = size(H_t);
[~,H_x_SIZE] = size(H_x);
H_t = repmat(H_t',1,H_x_SIZE);
H_x = repmat(H_x,H_t_SIZE,1);
H = H_t+H_x;
phantom_real = repmat(phantom,1500,1);

%%
out = MPI_sim_1D_1(H,f_s,phantom_real);
[~,L] = size(out);
% 截取一个周期并降采样
out = out(round(t*f_s+1):round(2*t*f_s));
out= downsample(out,10);
end