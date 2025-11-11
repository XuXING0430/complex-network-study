%% 复杂网络 SIR 模型仿真（MATLAB 可视化示例）
clc; clear; close all;

%% 网络参数
N = 50;           % 节点数
p = 0.07;         % ER 网络连边概率
steps = 30;       % 仿真步数

%% 构建随机网络（ER 网络）
G = rand(N) < p;
G = triu(G,1);    % 只保留上三角
G = G + G';       % 对称邻接矩阵

%% 初始化节点状态
state = repmat("S", N, 1); % 全部初始为易感者
state(randi(N)) = "I";     % 随机选择一个感染者

%% SIR 模型参数
beta = 0.3;   % 感染率
gamma = 0.1;  % 康复率

%% 节点布局（圆形）
theta = linspace(0, 2*pi, N+1);
pos = [cos(theta(1:N))', sin(theta(1:N))'];

%% 定义颜色映射
color_map = containers.Map(["S","I","R"], ...
                          {[0.7 0.7 0.7], ...  % S 灰色
                           [1 0.2 0.2], ...    % I 红色
                           [0.2 0.8 0.2]});    % R 绿色

%% 仿真循环
figure('Position',[100 100 600 600]);
for t = 1:steps
    new_state = state;
    
    % 状态更新
    for i = 1:N
        if state(i) == "I"
            neighbors = find(G(i,:));
            for j = neighbors
                if state(j) == "S" && rand < beta
                    new_state(j) = "I";
                end
            end
            if rand < gamma
                new_state(i) = "R";
            end
        end
    end
    state = new_state;
    
    % 绘制网络
    clf;
    gplot(G, pos, 'k-');  % 绘制连边
    hold on;
    
    % 每个节点的颜色
    node_colors = zeros(N,3);
    for k = 1:N
        node_colors(k,:) = color_map(state(k));
    end
    
    scatter(pos(:,1), pos(:,2), 100, node_colors, 'filled');
    title(['SIR 网络仿真, t = ' num2str(t)]);
    axis equal;
    axis off;
    drawnow;
end
