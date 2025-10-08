% 复杂网络指标计算主函数
% 作者：Three Purple
% 功能：从Excel读取邻接矩阵，计算网络拓扑指标

clc; clear;

%% 1. 读取邻接矩阵
filename = 'network.xlsx'; % Excel文件名
A = xlsread(filename);     % 读取邻接矩阵

% 检查矩阵有效性
if ~ismatrix(A) || size(A,1) ~= size(A,2)
    error('邻接矩阵必须是方阵');
end

N = size(A,1);             % 节点数
M = sum(A(:))/2;           % 无向网络边数（上三角求和更准确）
fprintf('节点数 N = %d\n', N);
fprintf('边数 M = %d\n', M);

%% 2. 网络密度
D = 2*M / (N*(N-1));
fprintf('网络密度 D = %.4f\n', D);

%% 3. 平均路径长度与直径
% 使用 graph 对象求最短路径
G = graph(A);
distMatrix = distances(G); 
L = mean(distMatrix(~isinf(distMatrix) & distMatrix>0), 'all'); % 平均最短路径
D_max = max(distMatrix(~isinf(distMatrix)), [], 'all');          % 直径
fprintf('平均路径长度 L = %.4f\n', L);
fprintf('网络直径 D_max = %.4f\n', D_max);

%% 4. 平均聚类系数
C_nodes = clustering_coef_bu(A);  % 调用辅助函数（见下方）
C = mean(C_nodes);
fprintf('平均聚类系数 C = %.4f\n', C);

%% 5. 节点中心性指标
k = sum(A,2);                           % 度
deg_centrality = k / (N-1);             % 度中心性
bet_centrality = centrality(G,'betweenness'); % 介数中心性
close_centrality = centrality(G,'closeness'); % 接近中心性
pagerank_score = centrality(G,'pagerank');    % PageRank

%% 6. 模块度与社区划分
% 使用Louvain算法（MATLAB R2021b+ 支持 communityDetection）
try
    comm = community_louvain(A);   % 需要有Brain Connectivity Toolbox (BCT)
    Q = modularity_louvain(A,comm);
    fprintf('模块度 Q = %.4f\n', Q);
catch
    warning('未检测到BCT工具箱，跳过模块度计算');
    comm = [];
    Q = NaN;
end

%% 7. 结果汇总输出
metrics = table((1:N)', k, deg_centrality, bet_centrality, close_centrality, pagerank_score, ...
    'VariableNames', {'Node','Degree','DegreeCentrality','Betweenness','Closeness','PageRank'});

disp('--- 节点级指标 ---');
disp(metrics);

fprintf('\n平均聚类系数 = %.4f\n', C);
fprintf('平均路径长度 = %.4f\n', L);
fprintf('网络直径 = %.4f\n', D_max);
fprintf('网络密度 = %.4f\n', D);
fprintf('模块度 Q = %.4f\n', Q);

%% 8. 可选：绘图
figure;
plot(G);
title('网络结构图');

function C = clustering_coef_bu(A)
% 计算无向无权网络的局部聚类系数
N = size(A,1);
C = zeros(N,1);
A = A - diag(diag(A)); % 移除自环
for i = 1:N
    neighbors = find(A(i,:));
    k = numel(neighbors);
    if k >= 2
        subA = A(neighbors, neighbors);
        E = sum(subA(:)) / 2;
        C(i) = 2*E / (k*(k-1));
    else
        C(i) = 0;
    end
end
end