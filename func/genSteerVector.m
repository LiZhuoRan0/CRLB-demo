function a = genSteerVector(theta, N, d, lambda)
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
a = exp(1j*2*pi*d*theta/lambda*n);

end