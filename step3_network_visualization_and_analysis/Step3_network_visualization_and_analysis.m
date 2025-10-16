% 复杂网络系列③：网络可视化与基础统计分析
% 作者: threepurple
% 日期: 2025-10-16

clc; clear;

%% 一、读取或构建网络
if isfile('edge_list.csv')
    T = readtable('edge_list.csv');
    G = graph(T.source, T.target);
    disp('成功导入边列表文件 edge_list.csv');
else
    disp('未检测到 edge_list.csv，生成小世界网络示例...');
    N = 100; K = 4; beta = 0.2;
    G = WattsStrogatz(N, K, beta);
end

%% 二、网络可视化
figure('Name','网络可视化');
p = plot(G, 'Layout', 'force');
title('网络可视化（Force Layout）');

% 根据节点度调整节点颜色与大小
deg = degree(G);
p.NodeCData = deg;
p.MarkerSize = 5 + deg;
colorbar;
title('网络可视化（节点颜色代表度）');

%% 三、基础统计分析
numNodes = numnodes(G);
numEdges = numedges(G);
fprintf('节点数: %d, 边数: %d\n', numNodes, numEdges);

% 平均度与分布
avgDeg = mean(deg);
fprintf('平均度: %.2f\n', avgDeg);
figure('Name','度分布');
histogram(deg);
xlabel('Degree'); ylabel('Frequency');
title('节点度分布');

%% 四、聚类系数与平均最短路径长度
C = mean(clusteringCoefficient(G));
L = mean(distances(G), 'all', 'omitnan');
fprintf('平均聚类系数: %.4f\n平均最短路径长度: %.4f\n', C, L);

%% 五、小世界特性分析
% 生成与G相同节点数和边数的随机网络
G_rand = randomReference(G);

C_rand = mean(clusteringCoefficient(G_rand));
L_rand = mean(distances(G_rand), 'all', 'omitnan');

SWI = (C / C_rand) / (L / L_rand);
fprintf('小世界指数 SWI = %.3f\n', SWI);
if SWI > 1
    disp('该网络表现出小世界特性');
else
    disp('该网络不具备显著小世界特性');
end

%% 六、无标度网络与幂律分析
deg = degree(G);
[counts, bins] = histcounts(deg, 'BinMethod', 'integers');
bins_center = bins(1:end-1) + diff(bins)/2;

figure('Name','幂律度分布');
loglog(bins_center, counts/sum(counts), 'o');
xlabel('k'); ylabel('P(k)');
title('度分布（对数坐标）');

% 幂律拟合（线性近似）
valid = bins_center > 0 & counts > 0;
logk = log(bins_center(valid));
logP = log(counts(valid)/sum(counts));
coeffs = polyfit(logk, logP, 1);
gamma = -coeffs(1);
fprintf('拟合幂律指数 γ ≈ %.3f\n', gamma);

hold on;
loglog(bins_center(valid), exp(polyval(coeffs, logk)), 'r--', 'LineWidth', 1.2);
legend('数据', sprintf('幂律拟合 (γ = %.2f)', gamma));
hold off;

%% 七、综合结论
disp('------ 网络特性总结 ------');
if SWI > 1
    disp('✅ 小世界特性：存在');
else
    disp('❌ 小世界特性：不显著');
end
if gamma > 2 && gamma < 3
    disp('✅ 无标度特性：幂律分布合理 (γ ∈ [2,3])');
else
    disp('❌ 无标度特性：未显著表现幂律分布');
end

%% 附加函数定义
function C = clusteringCoefficient(G)
    A = adjacency(G);
    n = size(A,1);
    C = zeros(n,1);
    for i=1:n
        ki = sum(A(i,:));
        if ki > 1
            subA = A(i,:) * A * A(:,i);
            C(i) = subA / (ki*(ki-1));
        end
    end
end

function G_rand = randomReference(G)
    % 构建具有相同节点数与边数的随机网络
    N = numnodes(G);
    E = numedges(G);
    edges = zeros(0,2); % 初始化为 0x2 矩阵，避免 ismember 报错
    while size(edges,1) < E
        s = randi(N);
        t = randi(N);
        if s ~= t && ~any(ismember(edges, [s t], 'rows')) && ~any(ismember(edges, [t s], 'rows'))
            edges = [edges; s t];
        end
    end
    G_rand = graph(edges(:,1), edges(:,2));
end
