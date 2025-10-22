%% step4_community_analysis.m
% ��������������ŵ���ͨ��4�� ���Žṹ����
% ���ߣ�three purple

clear; clc; close all;

%% 1. ����򹹽�����
if exist('adj_matrix.xlsx','file')
    A = readmatrix('adj_matrix.xlsx');
    % ȷ������Գơ����Ի�
    A = triu(A,1);
    A = A + A';
    A(A > 0) = 1;  % ��ֵ��������Ȩ�ɱ�����
    G = graph(A);
    disp('? �Ѵ� adj_matrix.xlsx �����ڽӾ���');
else
   N = 100; K = 4; beta = 0.2;
   G = WattsStrogatz(N, K, beta);
   A = adjacency(G);  % �õ�ϡ���ڽӾ���
    A = full(A); 
end
  
%% 2. Louvain ���ż�⣨����ʹ�� community_louvain��
if exist('community_louvain.m','file')
    [community, Q] = community_louvain(full(adjacency(G)));
else
    warning('δ�ҵ� community_louvain.m��ʹ������ genlouvain ʵ�֡�');
    [community, Q] = louvain_local(A);
end

fprintf('ģ��� Q = %.4f\n', full(Q));

%% 3. ���ӻ�����
figure;
p = plot(G, 'Layout', 'force');
title('�������Ż��ֽ��');
colors = lines(max(community));
for i = 1:max(community)
    highlight(p, find(community==i), 'NodeColor', colors(i,:));
end

%% 4. �������ͳ��
numCommunities = max(community);
fprintf('����������%d\n', numCommunities);
for i = 1:numCommunities
    fprintf('���� %d �ڵ�����%d\n', i, sum(community==i));
end

%% 5. �Žڵ����
bet = centrality(G, 'betweenness');
[~, idx] = max(bet);
fprintf('�Žڵ��ţ�%d�����������ԣ�%.3f\n', idx, bet(idx));

%% 6. �������ż����Ӿ���
commMatrix = zeros(numCommunities);
for i = 1:numCommunities
    for j = 1:numCommunities
        commMatrix(i,j) = sum(sum(A(community==i, community==j)));
    end
end
figure;
imagesc(commMatrix);
colorbar;
title('���ż�����ǿ�Ⱦ���');
xlabel('���ű��'); ylabel('���ű��');

%% 7. �Žڵ���ӻ�
figure;
p = plot(G, 'Layout', 'force');
highlight(p, idx, 'NodeColor', 'k', 'Marker', 's', 'MarkerSize', 8);
title('�Žڵ����');

%% ============================================================
% ��Ƕ Louvain ģ����Ż��������������ⲿ����ʱ���ã�
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