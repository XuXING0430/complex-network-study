% 复杂网络入门到精通7：网络鲁棒性基础 代码演示版本
clc; clear;

% 平均攻击 随机次数
repeat = 20;          
useApproxBet = true; 

X = 1:215;
% 调用主逻辑
results = core_attack_logic(repeat, useApproxBet);

NetEff_init = results.NetEff_init;
Max_Sub_init = results.Max_Sub_init;

%% 
% 图1: 最大连通子图规模 (归一化)
figure(1); hold on; box on;
plot([0 X], [1, results.Max_Sub_Random_aver/Max_Sub_init], '-b>', 'LineWidth',1.3,'MarkerSize',4);
plot([0 X], [1, results.Max_Sub_Degree/Max_Sub_init], '-gs', 'LineWidth',1.3,'MarkerSize',4);
plot([0 X], [1, results.Max_Sub_Bet/Max_Sub_init], '-ko', 'LineWidth',1.3,'MarkerSize',4);
legend('随机攻击','度攻击','介数攻击(近似)','Location','best','FontName','宋体');
xlabel('删除节点数量'); ylabel('最大连通子图 (归一)');
set(gca,'FontSize',14,'FontName','宋体');

% 图2: 全局效率 (归一化)
figure(2); hold on; box on;
plot([0 X], [1, results.NetEff_Random_aver/NetEff_init], '-b>', 'LineWidth',1.3,'MarkerSize',4);
plot([0 X], [1, results.NetEff_Degree/NetEff_init], '-gs', 'LineWidth',1.3,'MarkerSize',4);
plot([0 X], [1, results.NetEff_Bet/NetEff_init], '-ko', 'LineWidth',1.3,'MarkerSize',4);
legend('随机攻击','度攻击','介数攻击(近似)','Location','best','FontName','宋体');
xlabel('删除节点数量'); ylabel('全局效率 (归一)');
set(gca,'FontSize',14,'FontName','宋体');