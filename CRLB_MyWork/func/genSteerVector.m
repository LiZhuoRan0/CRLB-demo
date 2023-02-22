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
n = ((2*(0:N-1)-N+1)/2).';% ����ʦ���ɷ�ʽ���ǶȰ��м����Ԫ��
% n = 1:N;
a = zeros(N, length(theta));
for i = 1:length(theta)
%     a(:, i) = exp(1j*2*pi*1/2*sin(theta(i))*(n - 1)).';% ����ʸ��,û�й�һ��
%     a(:, i) = exp(1j*2*pi*d/lambda*theta*(n - 1)).'/sqrt(N);% ����ʸ��,��һ��
    a(:, i) = exp(1j*2*pi*d/lambda*theta*n).'/sqrt(N);% ����ʸ��,��һ��
%     a(:, i) = exp(1j*2*pi*1/2*theta*(n - 1)).'/sqrt(N);% ����ʸ��,��һ��
end

end