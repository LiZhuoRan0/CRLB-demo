%% 参数设定
clear
clc
addpath('.\func')
fprintf('==============================================\n');
rng(666);
%% 参数设置
N               =           256;                % 基站处天线
NRIS            =           256;                % RIS阵元数
NRF             =           4;
fc              =           100e9;              %100GHz
lambda_c        =           3e8/fc;
d               =           lambda_c/2;
B               =           10e9;               %10GHz
M               =           2048;                 %子载波数
angleOverSmp    =           2;
% S_eff           =           2;%RIS散射的面积
S_eff           =           (N*d)^2;
L_2             =           [5 1];
threshold       =           [0.85 0.95];

P_BS                           = 16;% 观测帧数
P_RIS                           = 32;

Nit             =           10;
% Nit             =           5;
L_max           =           20;%OMP最大迭代次数
%% -    BS，RIS，UE, scatter位置设定
% RIS位置
RIS.x                          = 20*sqrt(2);
RIS.y                          = 0;

angle_RIS                    = -pi/4;% BS ULA 朝向
RIS.x1                         = RIS.x - d*(NRIS-1)/2*abs(sin(angle_RIS));
RIS.y1                         = RIS.y + d*(NRIS-1)/2*abs(sin(angle_RIS));
RIS.xN                         = RIS.x + d*(NRIS-1)/2*abs(sin(angle_RIS));
RIS.yN                         = RIS.y - d*(NRIS-1)/2*abs(sin(angle_RIS));

% BS位置
BS.x                          = 0;% BS中心横纵坐标
BS.y                          = 0;
angle_BS                    = pi/4;% BS ULA 朝向
BS.x1                         = BS.x - d*(N-1)/2*abs(sin(angle_BS));% 以左边为第一阵元
BS.xN                         = BS.x + d*(N-1)/2*abs(sin(angle_BS));
BS.y1                         = BS.y - d*(N-1)/2*abs(sin(angle_BS));
BS.yN                         = BS.y + d*(N-1)/2*abs(sin(angle_BS));
% UE位置
UE.x                         = 20.1;
UE.y                         = -10.1;
% UE.x                         = RIS.x / 2;
% UE.y                         = -5;

r_UE2BS                     = sqrt((UE.x-BS.x).^2 +...
                                    (UE.y-BS.y).^2);                                
r_index_UE2BS               = r_UE2BS*B/3e8;
r_UE2RIS                    = sqrt((UE.x-RIS.x).^2 +...
                                    (UE.y-RIS.y).^2); 
r_index_UE2RIS               = r_UE2RIS*B/3e8; 
r_BS2RIS                    = sqrt((BS.x - RIS.x).^2 + (BS.y - RIS.y).^2);
theta_UE2BS                 = atan((UE.y-BS.y)/(UE.x-BS.x));         
theta_UE2RIS                = atan((UE.y-RIS.y)/(UE.x-RIS.x));

%% -    发射功率，信噪比,未经由RIS
P_noise = 10^(-17.4)*(60e3)*(M);
P_noise_dBm = 10*log10(P_noise);

P_t_dBm = 0:4:20;%dBmW
% P_t_dBm = 0;
CRLB = zeros(length(P_t_dBm), 1);
CRLB_r = zeros(length(P_t_dBm), 1);
NMSE = zeros(length(P_t_dBm), length(L_2));
Times = zeros(length(P_t_dBm), length(L_2));
G_r = 1;
G_t = 1;
P_r_dBm = P_t_dBm + 10*log10(G_r*G_t*lambda_c^2/(4*pi)^2/(r_UE2BS)^2);

SNRdBs = P_r_dBm - P_noise_dBm;
%% -    经由RIS
G_r_RIS = N;
P_r_RIS_dBm = P_t_dBm + 10*log10(G_r_RIS*G_t*S_eff*lambda_c^2/(4*pi)^3/(r_UE2RIS)^2/(r_BS2RIS)^2);
SNRdBs_RIS = P_r_RIS_dBm - P_noise_dBm;
%% -    -   散射体
num_scatter_UE2BS = 5;
num_scatter_UE2RIS= 5;
num_scatter_com = 0;
PathInCluster = 6;
SizeOfScatter = 1;
% 公共散射体
Scatter_com.x = 10+(20-10)*rand(1, num_scatter_com);
Scatter_com.y = -10+(-1-(-10))*rand(1, num_scatter_com);
% UE2BS散射体
if num_scatter_com ~= 0
    if num_scatter_UE2BS ~= 0
        Scatter_UE2BS.x = [Scatter_com.x  5+(20-5)*rand(1, num_scatter_UE2BS)];   
        Scatter_UE2BS.y = [Scatter_com.y  -25+(-5-(-25))*rand(1, num_scatter_UE2BS)];
    else
        Scatter_UE2BS.x = [Scatter_com.x];   
        Scatter_UE2BS.y = [Scatter_com.y];
    end
else
    if num_scatter_UE2BS ~= 0
        Scatter_UE2BS.x = [5+(20-5)*rand(1, num_scatter_UE2BS)];   
        Scatter_UE2BS.y = [-25+(-5-(-25))*rand(1, num_scatter_UE2BS)];
    else
        Scatter_UE2BS.x = [];   
        Scatter_UE2BS.y = [];
    end
end
if (num_scatter_UE2BS+num_scatter_com) ~= 0  
    Scatter_UE2BS_tmp.x = Scatter_UE2BS.x;
    Scatter_UE2BS_tmp.y = Scatter_UE2BS.y;
    Scatter_UE2BS.x = [];
    Scatter_UE2BS.y = [];
    % 在生成的Scatter_UE2BS上加簇
    for i = 1:PathInCluster
        Scatter_UE2BS.x = [Scatter_UE2BS.x, Scatter_UE2BS_tmp.x-SizeOfScatter/2+SizeOfScatter/PathInCluster*i];
        Scatter_UE2BS.y = [Scatter_UE2BS.y, Scatter_UE2BS_tmp.y];
    end
end
Scatter_UE2BS.theta = atan((Scatter_UE2BS.y-BS.y)./(Scatter_UE2BS.x-BS.x));
Scatter_UE2BS.r = distance(repmat(UE.x, 1, length(Scatter_UE2BS.x))...
                , repmat(UE.y, 1, length(Scatter_UE2BS.y))...
                , Scatter_UE2BS.x, Scatter_UE2BS.y)...
                + distance(repmat(BS.x, 1, length(Scatter_UE2BS.x))...
                , repmat(BS.y, 1, length(Scatter_UE2BS.y))...
                , Scatter_UE2BS.x, Scatter_UE2BS.y);
% Scatter_UE2BS.alpha = 1./(Scatter_UE2BS.r);% 经过散射体幅度下降6dB
Scatter_UE2BS.alpha = 1./distance(repmat(UE.x, 1, length(Scatter_UE2BS.x))...
                , repmat(UE.y, 1, length(Scatter_UE2BS.y))...
                , Scatter_UE2BS.x, Scatter_UE2BS.y)...
                .*1./distance(repmat(BS.x, 1, length(Scatter_UE2BS.x))...
                , repmat(BS.y, 1, length(Scatter_UE2BS.y))...
                , Scatter_UE2BS.x, Scatter_UE2BS.y)...
                ./sqrt(PathInCluster);

% UE2RIS散射体
if num_scatter_com ~= 0
    if num_scatter_UE2RIS ~= 0        
        Scatter_UE2RIS.x = [Scatter_com.x  5+(20-5)*rand(1, num_scatter_UE2RIS)];   
        Scatter_UE2RIS.y = [Scatter_com.y  -25+(-5-(-25))*rand(1, num_scatter_UE2RIS)];
    else
        Scatter_UE2RIS.x = [Scatter_com.x];   
        Scatter_UE2RIS.y = [Scatter_com.y];
    end
else
    if num_scatter_UE2RIS ~= 0
        Scatter_UE2RIS.x = [5+(20-5)*rand(1, num_scatter_UE2RIS)];   
        Scatter_UE2RIS.y = [-25+(-5-(-25))*rand(1, num_scatter_UE2RIS)];
    else
        Scatter_UE2RIS.x = [];   
        Scatter_UE2RIS.y = [];
    end
end
if (num_scatter_UE2RIS+num_scatter_com) ~= 0  
    Scatter_UE2RIS_tmp.x = Scatter_UE2RIS.x;
    Scatter_UE2RIS_tmp.y = Scatter_UE2RIS.y;
    Scatter_UE2RIS.x = [];
    Scatter_UE2RIS.y = [];
    % 在生成的Scatter_UE2BS上加簇
    for i = 1:PathInCluster
        Scatter_UE2RIS.x = [Scatter_UE2RIS.x, Scatter_UE2RIS_tmp.x-SizeOfScatter/2+SizeOfScatter/PathInCluster*i];
        Scatter_UE2RIS.y = [Scatter_UE2RIS.y, Scatter_UE2RIS_tmp.y];
    end
end
Scatter_UE2RIS.theta = atan((Scatter_UE2RIS.y-RIS.y)./(Scatter_UE2RIS.x-RIS.x));
Scatter_UE2RIS.r = distance(repmat(UE.x, 1, length(Scatter_UE2RIS.x))...
                , repmat(UE.y, 1, length(Scatter_UE2RIS.y))...
                , Scatter_UE2RIS.x, Scatter_UE2RIS.y)...                
                + distance(repmat(RIS.x, 1, length(Scatter_UE2RIS.x))...
                , repmat(RIS.y, 1, length(Scatter_UE2RIS.y))...
                , Scatter_UE2RIS.x, Scatter_UE2RIS.y);
Scatter_UE2RIS.alpha = 1./(distance(repmat(UE.x, 1, length(Scatter_UE2RIS.x))...
                , repmat(UE.y, 1, length(Scatter_UE2RIS.y))...
                , Scatter_UE2RIS.x, Scatter_UE2RIS.y))...
                ...
                .*1./(distance(repmat(RIS.x, 1, length(Scatter_UE2RIS.x))...
                , repmat(RIS.y, 1, length(Scatter_UE2RIS.y))...
                , Scatter_UE2RIS.x, Scatter_UE2RIS.y))...
                ./sqrt(PathInCluster);% 经过散射体幅度下降6dB
%% 参数输出
fprintf('P_BS = %d\n', P_BS);
fprintf('P_RIS = %d\n', P_RIS);
fprintf('num_scatter_UE2BS = %d\n', num_scatter_UE2BS);
fprintf('num_scatter_UE2RIS = %d\n', num_scatter_UE2RIS);
fprintf('num_scatter_com = %d\n', num_scatter_com);
fprintf('真实的UE2BS距离 = %10.6f\n', r_UE2BS);                            
fprintf('真实的UE2RIS距离 = %10.6f\n', r_UE2RIS);
fprintf('真实的UE2BS角度（弧度） = %10.8f\n', theta_UE2BS); 
fprintf('真实的UE2BS角度（相对BS的法线，无单位） = %10.8f\n', sin(pi/2-angle_BS+theta_UE2BS)); 
fprintf('真实的UE2RIS角度（弧度） = %10.8f\n', theta_UE2RIS); 
fprintf('真实的UE2RIS角度（相对RIS的法线，无单位） = %10.8f\n', sin(pi/2+angle_RIS-theta_UE2RIS)); 
%% Y和A生成
H_UE2RIS_absolute = NearFieldH(NRIS, sin(pi/2+angle_RIS-[theta_UE2RIS, Scatter_UE2RIS.theta])...            
            , [r_UE2RIS, Scatter_UE2RIS.r], [1/(r_UE2RIS), Scatter_UE2RIS.alpha], fc, B, M, 'absolute');
H_BR = zeros(N, NRIS, M);
for i_M = 1:M
   H_BR(:, :, i_M) = genRIS2BS_H(RIS.x1, RIS.y1,...
                        RIS.xN, RIS.yN,...
                        BS.x1, BS.y1,...
                        BS.xN, BS.yN,...
                        NRIS, N, fc-B/2+(B/M)*(i_M-1))./r_BS2RIS; %fc-B/2+(B/M)*(i_M-1)
end

H_UE2BS_absolute = NearFieldH(N, sin(pi/2-angle_BS+[theta_UE2BS, Scatter_UE2BS.theta])...
            , [r_UE2BS, Scatter_UE2BS.r], [1/(r_UE2BS) Scatter_UE2BS.alpha], fc, B, M, 'absolute');
HB  = H(1/r_UE2BS, fc, r_UE2BS, N, sin(pi/2-angle_BS+theta_UE2BS), M, B);
for i_P_t_dBm = 1:length(P_t_dBm)
    fprintf('==============================================\n');        
    fprintf('P_t_dBm = %3d\n',P_t_dBm(i_P_t_dBm));
    for i_it = 1:Nit 
        if mod(i_it, 10) == 0
            fprintf('i_it/Nit = %3d/%3d\n', i_it, Nit);
        end
%% -    根据参数生成Y，A，H
%% -    -   UE2BS        
        Y_absolute_all = zeros(NRF*P_BS, M);
        A = exp(1j*2*pi*rand(NRF*P_BS, N));
        A(1,:) = zeros(1, N);
        A(1, ceil(N/2)) = 1;
        A(1, ceil((N+1)/2)) = 1;
        
        for i_P = 1:P_BS
            for i_NRF = 1:NRF
                Y_absolute_all(i_NRF + (i_P-1)*NRF, :) = awgn(A(i_NRF + (i_P-1)*NRF, :)*H_UE2BS_absolute, SNRdBs(i_P_t_dBm), 'measured');
            end
        end
        %% CRLB                
        %% -    求Fisher矩阵
        sigma2  = norm(Y_absolute_all, 'fro')^2 /...
                (size(Y_absolute_all, 1)*size(Y_absolute_all, 2)) / ...
                (10^(SNRdBs(i_P_t_dBm)/10)+1);
        FB  = FisherB(Y_absolute_all, A, HB, sigma2, H_UE2BS_absolute);
        tem = sqrt(inv(FB));
        tem_t = tem(2,2);
        tem_r = tem(3,3);
        CRLB(i_P_t_dBm) = CRLB(i_P_t_dBm) + tem_t / (1-sin(pi/2-angle_BS+theta_UE2BS)^2);
        CRLB_r(i_P_t_dBm) = CRLB_r(i_P_t_dBm) + tem_r;
    end
    CRLB(i_P_t_dBm) = CRLB(i_P_t_dBm)/Nit;
    CRLB_r(i_P_t_dBm) = CRLB_r(i_P_t_dBm)/Nit;
    fprintf('UE2BS 角度(弧度) CRLB = %5.4f 1e-6\n', CRLB(i_P_t_dBm)*1e6);
    fprintf('UE2BS 距离(米)   CRLB = %5.4f 1e-4\n', CRLB_r(i_P_t_dBm)*1e4);
end 
figure; hold on; box on; grid on;
plot(P_t_dBm, real(CRLB), 'ro-');
% plot(P_t_dBm, real(CRLB_r), 'r<--');
set(gca, 'YScale', 'log')
set(gca,'YLim',[1e-8 1e-3]);
legend('CRLB_{\vartheta}')
save('./Data/CDL_Near_2048', 'P_t_dBm', 'CRLB')