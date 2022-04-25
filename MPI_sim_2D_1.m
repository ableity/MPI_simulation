function [out_x,out_y] = MPI_sim_2D_1(H_x,H_y,f_s,phantom)
%% 二维无弛豫MPI仿真函数（由MPI_sim_2D_2调用）
%% 
% 李蕾 2022年04月19日
% 二维MPI仿真程序，修改自一维
% 建立一个矩阵，一行数据是在相同时间不同位置下的H
% 一列数据是在相同位置下不同时间的t
% 总结起来就是行代表位置，列代表时间
% 那么，无论是几维的数据，都可以通过拉成一列数据来模拟，类似于系统矩阵
%% 第一部分，输入参数

%磁粒子直径，单位为nm
D = 20;
%接收灵敏度，单位为mT/A
S = 9;
%溶液体积，单位为ml
VT = 0.2;
%溶液铁浓度，单位mg/mL
n_Fe = 25;
%弛豫时间 （微秒），本程序没有用到
relaxation_time = 1;
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
%弛豫时间单位转化，本程序没有用到
relaxation_time = relaxation_time*1e-6;

%参数集合，方便后续输入
parameter_M_H.c = c;
parameter_M_H.m = m;
parameter_M_H.B = B;
parameter_M_H.u0 = u0;
parameter_M_H.VT = VT;
parameter_M_H.S = S;

%% 第四部分，外部参数，已经改为在第二个外部函数中编辑
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
% 粒子磁矩，H_x和H_y是行为位置，列为时间的矩阵，矩阵中存储外加磁场大小
% 首先将两个方向合成一个方向
H = (H_x.^2+H_y.^2).^0.5;
%按照仿体稀疏性加速，实测这儿只能加速2s左右
%取出仿体的一行（实际每行都是一样的）
phantom_single = phantom(1,:);
%找出为0的位置
phantom_num_zero = find(phantom_single==0);
%找出非0的位置
phantom_num_not_zero = find(phantom_single~=0);
%准备裁剪的仿体矩阵
phantom_not_zero = phantom;
%准备裁剪的外磁场矩阵
H_not_zero = H;
%将仿体0的位置置空，即不需要计算
phantom_not_zero(:,phantom_num_zero)=[];
%将外磁场为0的位置置空，即不需要计算
H_not_zero(:,phantom_num_zero)=[];
%准备一个矩阵准备记录M_p
M_p_total = zeros(size(H));
% MPI仿真，按照郎之万模型

M_p_not_zero = M_H(parameter_M_H,H_not_zero,phantom_not_zero);
%补足之前的零位置
M_p_total(:,phantom_num_zero) = 0;
%补足非0位置

M_p_total(:,phantom_num_not_zero) = M_p_not_zero;
M_p = M_p_total;
% 将M分解为X和Y两个方向，因为接收线圈只能接收一个方向
M_p_x = M_p.*H_x./H;
M_p_y = M_p.*H_y./H;
% 直接按照上文计算的外磁场H比例分解方向，但是会出现0/0而等于NaN的情况
% 按照程序逻辑NaN处的值均应该为0，所以找到这些值，将其置零
M_p_x(isnan(M_p_x))=0;
M_p_y(isnan(M_p_y))=0;

% 将空间维度累加，因为接收线圈并不能区分不同位置的信号
M_p_x = sum(M_p_x');
M_p_y = sum(M_p_y');

%磁矩求导，接收的信号，求导后即输出
M_p_d_x = [M_p_x(2:end),M_p_x(end)]-M_p_x;
M_p_d_x = M_p_d_x/(1/f_s);

M_p_d_y = [M_p_y(2:end),M_p_y(end)]-M_p_y;
M_p_d_y = M_p_d_y/(1/f_s);

%（3）后处理，与预处理对应
% [~,L] = size(M_p_d);
% % 截取一个周期并降采样
% M_p_d_3 = M_p_d;
% M_p_d_3 = M_p_d_3(round(t*f_s+1):round(2*t*f_s));
% M_p_d= downsample(M_p_d_3,10);
out_y = M_p_d_y;
out_x = M_p_d_x;
%% 函数
function out = M_H(parameter,in,phantom)
    % 加参数的郎之万函数
    out = parameter.c*parameter.m*langevin(parameter.B*in/parameter.u0);
    out = phantom.*out;
    out = parameter.u0*parameter.S*parameter.VT*out;
end

function out=langevin(in)
    % 未加参数的郎之万函数，0.000001用来防止分母为0
    in = in+0.000001;
    out = coth(in)-1./in;
end


end
