function out = MPI_sim_1D_1(H,f_s,phantom)
%% 一维弛豫MPI仿真函数
%%
% 李蕾 2022年04月19日
% 二维MPI仿真程序，修改自一维
% 一维MPI的磁场输入为
% 建立一个矩阵，一行数据是在相同时间不同位置下的H
% 一列数据是在相同位置下不同时间的t
% 总结起来就是行代表位置，列代表时间
% 那么，无论是几维的数据，都可以通过拉成一列数据来模拟，类似于系统矩阵
%% 第一部分，输入参数
%路径

%磁粒子直径，单位为nm
D = 20;
%接收灵敏度，单位为mT/A
S = 9;
%溶液体积，单位为ml
VT = 0.2;
%溶液铁浓度，单位mg/mL
n_Fe = 25;

%弛豫时间 （微秒）
relaxation_time = 0;

%% 第二部分，默认参数
%真空磁导率
u0=4*pi*1e-7;
%磁粒子饱和磁化强度(A/m为单位)
MS=446000;
%四氧化三铁密度(Kg/m^3)
p_Fe3O4 = 5200; 
%四氧化三铁摩尔质量(Kg/mol)
m_Fe3O4 = 0.232;
%玻尔兹曼常数
kB=1.38e-23;
%温度
Tp = 300;

%% 第三部分，计算参数，根据前面的参数对某些参数进行推导
%磁粒子直径转化单位为m
D = D*1e-9;
%磁粒子体积
V=1/6*pi*D^3;
%每个磁粒子核质量(Kg)
m_core = p_Fe3O4*V; 
%每个核中铁离子摩尔量(mol)
N_Fe = m_core/m_Fe3O4*3;
%磁矩
m=V*MS;
%计算出β作为系数
B=u0*m/(kB*Tp);
%计算溶液中四氧化三铁浓度（/L）
n_Fe = n_Fe/56*1e3;
c = n_Fe/N_Fe;
%灵敏度转化mT->A/m
S = S*1e-3/u0;
%溶液体积单位转化mL
VT = VT*1e-6;
%弛豫时间单位转化
relaxation_time = relaxation_time*1e-6;

%参数集合
parameter_M_H.c = c;
parameter_M_H.m = m;
parameter_M_H.B = B;
parameter_M_H.u0 = u0;
parameter_M_H.VT = VT;
parameter_M_H.S = S;

%% 第四部分，外部参数
% 本部分为外场所有参数，此部分可以单独写脚本当外置输入，所以分开写
%外加激励磁场峰值，单位mT
% H_peak = 12.5;
% %外加梯度磁场梯度,单位T/m
% H_gradient = 2.5;
% %峰值单位转化
% H_peak = H_peak*1e-3;
%采样频率
% f_s_raw = 1250000;
% %采样时间
% t = 1/25000;
% %激励频率
% f = 25000;
% %FOV细分数量
% num_of_pixel = 1000; 
% %FOV,一个长度，单位为m
% FOV = H_peak/H_gradient*2;
% FOV = FOV*1.2;
% %像素长度
% pixel = FOV/num_of_pixel;
% 
% % 增加采样率与采样点数
% % 时间序列（按照采样率的10倍，后面通过降采样回到真实采样率）
% % 时间用了3倍，后面可以取中间的一段
% 
% % 激励磁场
% f_s = 10*f_s_raw;
% t_real = 1/f_s:1/f_s:t*3;
% % 外加激励磁场随时间的变化
% H_t = H_peak*sin(2*pi*f*t_real);
% % 梯度磁场
% H_x = -0.5*H_gradient*FOV+pixel*H_gradient:pixel*H_gradient:0.5*H_gradient*FOV;
% 
% % 仿体数据
% phantom = zeros(1,1000);
% phantom(500)=1;
% phantom(1000)=1;
% 
% 
% % 建立一个矩阵，一行数据是在相同时间不同位置下的H
% % 一列数据是在相同位置下不同时间的t
% [~,H_t_SIZE] = size(H_t);
% [~,H_x_SIZE] = size(H_x);
% H_t = repmat(H_t',1,H_x_SIZE);
% H_x = repmat(H_x,H_t_SIZE,1);
% H = H_t+H_x;
% phantom_real = repmat(phantom,1500,1);
% H = H.*phantom_real;

% 下方在该部分需要用到的参数有
% H 存储磁场与粒子相对浓度信息
% f_s 求导需要时间
% t存储采样时间

%% （2）仿真
% 粒子磁矩（若有磁粒子下的单位磁矩）
M_p = M_H(parameter_M_H,H,phantom);
% 空间位置求和
M_p = sum(M_p');

% 弛豫
if relaxation_time~=0
    sample_point = floor(4*relaxation_time*f_s);
    t = (1:sample_point)/f_s;
    r = exp(-t/relaxation_time);
    r = r./sum(r);
    %策略1 先弛豫卷积再求导
     M_p_r = conv(M_p,r);
     M_p_r = M_p_r(1:size(M_p,2));
else
    M_p_r=M_p;
end

%磁矩求导，接收的信号
M_p_d_r = [M_p_r(2:end),M_p_r(end)]-M_p_r;
M_p_d_r = M_p_d_r/(1/f_s);



%（3）后处理，与预处理对应
% [~,L] = size(M_p_d);
% % 截取一个周期并降采样
% M_p_d_3 = M_p_d;
% M_p_d_3 = M_p_d_3(round(t*f_s+1):round(2*t*f_s));
% M_p_d= downsample(M_p_d_3,10);
out = M_p_d_r;

%% 函数
function out = M_H(parameter,in,phantom)
    out = parameter.c*parameter.m*langevin(parameter.B*in/parameter.u0);
    out = phantom.*out;
    
    out = parameter.u0*parameter.S*parameter.VT*out;
end

function out=langevin(in)
in = in+0.000001;
    out = coth(in)-1./in;
end
function out = fft_recode(in)
%将FFT的常规操作集合了起来

[a,b] = size(in);
if a==1
    L = b;
else
    L = a;
end
Y = fft(in,L);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
out = P1/(L/2);

end

function out = find_begin_point(in,L)
    % 找起始相位点
    % L为输出长度
    % 找到最大最最小值的平均值，作为起点
    average = (min(min(in))+max(max(in)))/2;
    % 寻找方法为，找到离该平均值最近的点
    temp = abs(in-average);
    % 负值寻峰
    begin = findpeaks(-temp);
    % 根据值寻找位置
    begin = find(begin(1) == -temp);
    %可能有多个点，取第一个
    begin = begin(1);
    %输出
    out = in(begin:begin+L-1);
end





end
