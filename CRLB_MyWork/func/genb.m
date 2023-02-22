function [b] = genb(theta, r, N, f, d)
%% ���ɽ�������ʸ��

% Inputs��
%       theta & r=          =       ��ʽ7�������
%       N                   =       ������
%       fc                  =       �ز�Ƶ��
%       d                   =       ������Ԫ���

% Outputs��
%       b                   =       ��������ʸ��
%%
% b = zeros(N, 1);
kc = 2*pi*f/3e8;
delta_n = ((2*(0:N-1)-N+1)/2).';% ����ʦ���ɷ�ʽ���ǶȰ��м����Ԫ��
% delta_n = (0:N-1).';% ��ʹ�õķ�ʽ���ǶȰ���һ����Ԫ��
% d = lambda/2;
r_l_n = sqrt(r^2 + delta_n.^2.*d.^2 - 2*r*theta*delta_n*d);

b = exp(-1j*kc*(r_l_n-r))/sqrt(N);
end