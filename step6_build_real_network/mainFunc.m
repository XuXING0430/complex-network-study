clc;clear;

% excel name
filename = 'test.xls';

% 输入 filename为文件名
% 此函数通过Excel文件 生成邻接矩阵A，邻接表edge_list，节点编号对应关系indexs
[A, edge_list, indexs] = buildNet(filename);