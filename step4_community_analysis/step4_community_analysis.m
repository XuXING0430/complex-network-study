%% step4_community_analysis.m
% 复杂网络仿真入门到精通（4） 社团结构分析
% 作者：three purple

clear; clc; close all;

%% 1. 载入或构建网络
if exist('adj_matrix.xlsx','file')
    A = readmatrix('adj_matrix.xlsx');
    % 确保矩阵对称、无自环
    A = triu(A,1);
    A = A + A';
    A(A > 0) = 1;  % 二值化（若含权可保留）
    G = graph(A);
    disp('? 已从 adj_matrix.xlsx 载入邻接矩阵');
else
   N = 100; K = 4; beta = 0.2;
   G = WattsStrogatz(N, K, beta);
   A = adjacency(G);  % 得到稀疏邻接矩阵
    A = full(A); 
end
  
%% 2. Louvain 社团检测（优先使用 community_louvain）
if exist('community_louvain.m','file')
    [community, Q] = community_louvain(full(adjacency(G)));
else
    warning('未找到 community_louvain.m，使用内置 genlouvain 实现。');
    [community, Q] = louvain_local(A);
end

fprintf('模块度 Q = %.4f\n', full(Q));

%% 3. 可视化社团
figure;
p = plot(G, 'Layout', 'force');
title('网络社团划分结果');
colors = lines(max(community));
for i = 1:max(community)
    highlight(p, find(community==i), 'NodeColor', colors(i,:));
end

%% 4. 输出社团统计
numCommunities = max(community);
fprintf('社团数量：%d\n', numCommunities);
for i = 1:numCommunities
    fprintf('社团 %d 节点数：%d\n', i, sum(community==i));
end

%% 5. 桥节点分析
bet = centrality(G, 'betweenness');
[~, idx] = max(bet);
fprintf('桥节点编号：%d，介数中心性：%.3f\n', idx, bet(idx));

%% 6. 绘制社团间连接矩阵
commMatrix = zeros(numCommunities);
for i = 1:numCommunities
    for j = 1:numCommunities
        commMatrix(i,j) = sum(sum(A(community==i, community==j)));
    end
end
figure;
imagesc(commMatrix);
colorbar;
title('社团间连接强度矩阵');
xlabel('社团编号'); ylabel('社团编号');

%% 7. 桥节点可视化
figure;
p = plot(G, 'Layout', 'force');
highlight(p, idx, 'NodeColor', 'k', 'Marker', 's', 'MarkerSize', 8);
title('桥节点高亮');

%% ============================================================
% 内嵌 Louvain 模块度优化函数（仅当无外部函数时启用）
% ============================================================
function [S, Q] = louvain_local(A)
    gamma = 1;
    k = full(sum(A));
    twom = sum(k);
    B = full(A - gamma * (k' * k) / twom);
    [S, Q] = genlouvain(B);
    Q = Q / twom;
    Q = full(Q);
end