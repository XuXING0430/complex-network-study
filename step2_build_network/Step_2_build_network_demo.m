%% build_network_demo.m
% 复杂网络构建演示脚本
% 作者: threepurple
% 日期: 2025-10-08

clc; clear;

%% 1. ER 随机网络
N = 10; p = 0.3;
A_er = rand(N) < p;
A_er = triu(A_er,1); A_er = A_er + A_er';
G_er = graph(A_er);
figure;
plot(G_er);
title('Erdős–Rényi 随机网络');

%% 2. 小世界网络
N = 20; K = 4; beta = 0.2;
G_ws = WattsStrogatz(N, K, beta);
A_ws = adjacency(G_ws);
figure;
plot(G_ws);
title('Watts–Strogatz 小世界网络');

%% 3. 无标度网络
N = 30; m0 = 3; m = 2;
A_ba = zeros(N);
A_ba(1:m0,1:m0) = 1 - eye(m0);
for i = (m0+1):N
    degrees = sum(A_ba(1:i-1,1:i-1));
    P = degrees / sum(degrees);
    targets = randsample(i-1, m, true, P);
    A_ba(i,targets) = 1; A_ba(targets,i) = 1;
end
G_ba = graph(A_ba);
figure;
plot(G_ba);
title('Barabási–Albert 无标度网络');

%% 4. 实体网络导入
% 假设存在 edge_list.csv 文件
if isfile('edge_list.csv')
    T = readtable('edge_list.csv');
    G_real = graph(T.source, T.target);
    figure; plot(G_real); title('实体网络 - 边列表导入');
end

disp('✅ 所有网络构建示例已完成');

