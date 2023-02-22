function [h] = gen1toArray(xp, yp,...
                           x1, y1,...
                           xN, yN,...
                           N, fc)
%% ����1�㵽һ�����е��ŵ�                       
% Inputs:
%       xp,yp�ǵ�Ŀ���λ��
%       x1, y1�������ͷ����
%       xN, yN�������β����
%       N:�������Ԫ��
% Outputs:
%       h                       =       ��Ŀ�굽������ŵ���N��1��
%%
lambda = 3e8/fc;
kc = 2*pi/lambda;
x = linspace(x1, xN, N).';% ���������x����
y = linspace(y1, yN, N).';% ���������y����
h = exp(-1j*kc* sqrt((x-xp).^2 + (y-yp).^2));

% [b] = genb(theta, r, N, fc);


end