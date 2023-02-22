function H = genRIS2BS_H(xR1, yR1,...
                        xRN, yRN,...
                        xB1, yB1,...
                        xBN, yBN,...
                        NRIS, NBS, fc)
%% ����RIS��BS���ŵ�                    
% Inputs:
%       RIS�ͻ�վ��ͷ��β������ֵ
%       RIS��BS��������
%       fc                              =       ��Ƶ
% Outputs��
%       RIS 2 BS ���ŵ�
%%
H = zeros(NBS, NRIS);
% ���������Ԫ�����������
xB = linspace(xB1, xBN, NBS);
yB = linspace(yB1, yBN, NBS);
% xR = linspace(xR1, xRN, NRIS);
% yR = linspace(yR1, yRN, NRIS);
for i_NBS = 1:NBS
    H(i_NBS, :) = gen1toArray(xB(i_NBS), yB(i_NBS),...
                   xR1, yR1,...
                   xRN, yRN,...
                   NRIS, fc).';
end

end