function a = genPartialSteerVector(theta, N, d, lambda, flag)
%%  ���ɵ���ʸ��,theta������Ƕ�

% inputs
%       theta                   =       ����Ƕ�
%       N                       =       ������Ԫ��
%       d                       =       ������Ԫ���
%       lambda                  =       ����

% outputs
%       a                       =       ���ɵĵ���ʸ����N��1��
%%
% n = ((2*(0:N-1)-N+1)/2).';% ����ʦ���ɷ�ʽ���ǶȰ��м����Ԫ��
n = (0:(N-1)).';
if flag == 1    
    a = (1j*pi.*n*cos(asin(theta))).*exp(1j*2*pi*d*theta/lambda*n);
else
    a = -(pi.*n.*cos(asin(theta))).^2.*exp(1j*2*pi*d*theta/lambda*n) -...
        (1j*pi.*n*theta).*exp(1j*2*pi*d*theta/lambda*n);
end
end