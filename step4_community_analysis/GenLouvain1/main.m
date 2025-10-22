% 社区算法 genlouvain

% 加载数据
load('adj_matrix.mat')

gamma = 1;
k = full(sum(A));
twom = sum(k);
B = full(A - gamma*k'*k/twom);
% 调用函数
[S,Q] = genlouvain(B);
% 模块度
Q = Q/twom;
% 查看有多少个分组
S_number = unique(S);