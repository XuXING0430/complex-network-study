% �����㷨 genlouvain

% ��������
load('adj_matrix.mat')

gamma = 1;
k = full(sum(A));
twom = sum(k);
B = full(A - gamma*k'*k/twom);
% ���ú���
[S,Q] = genlouvain(B);
% ģ���
Q = Q/twom;
% �鿴�ж��ٸ�����
S_number = unique(S);