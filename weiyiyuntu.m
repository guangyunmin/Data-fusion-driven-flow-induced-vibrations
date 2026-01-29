% ================== 零插值云图（修复少一根 + 红色边界线 + 字体修复） ==================
clear; clc; close all

%% 参数：你的数据文件名（空格或Tab分隔，无表头）
fname = 'data.txt';   % <- 修改为你的txt文件名

%% 读取与预处理
A = readmatrix(fname);                       % 自动识别科学计数法
if isempty(A) || size(A,2) < 2
    error('数据文件至少需要1列时间 + 1列位移。');
end

t_raw = A(:,1);                              % 时间（可不均匀）
U_raw = A(:,2:end);                          % 位移 (Ntime × Nrods)

% 按时间排序（若原始数据非严格递增）
[ t, idx ] = sort(t_raw);
U = U_raw(idx,:);

%% 组装网格（不插值）
Z = U.';                                     % 行=棒号，列=时间
[nRods, nTime] = size(Z);
if nRods == 0, error('没有检测到任何棒的位移列。'); end

% --- 关键：构造“边界坐标”并在 CData 末尾复制一行/一列，避免少一根 ---
% 时间边界（长度 = nTime+1），保持不均匀间隔的真实比例
if nTime > 1
    dt = diff(t);
    t_edges = [t(1)-dt(1)/2; (t(1:end-1)+t(2:end))/2; t(end)+dt(end)/2];
else
    t_edges = [t-0.5; t+0.5];                 % 仅一个时间点的兜底边界
end

% 棒号边界（长度 = nRods+1）
y_edges = (0.5 : 1 : nRods+0.5).';

% 将 Z 在末尾复制一行、一列
Ce = Z;
Ce(end+1,:) = Z(end,:);
Ce(:,end+1) = Ce(:,end);

% 构造与 Ce 同尺寸的坐标网格
[Xe, Ye] = meshgrid(t_edges, y_edges);

%% 作图
figure('Color','w');
h = pcolor(Xe, Ye, Ce);
shading flat
set(h, 'EdgeColor', 'none');

% === 在每根棒之间加边界线（红色） ===
hold on;
for k = 1:nRods+1
    plot([t_edges(1), t_edges(end)], [k-0.5, k-0.5], 'k-', 'LineWidth', 1.2); % 若要黑色：改 'r-' 为 'k-'
end
hold off;

%% 坐标轴与刻度
set(gca, 'YDir', 'normal', 'YLim', [0.5, nRods+0.5], 'YTick', 1:nRods);

% X 轴按 0.1 s 固定间隔
t_min = 0;         % 向下取整到 0.1 的倍数
t_max = 0.5;          % 向上取整到 0.1 的倍数
xticks = t_min:0.1:t_max;
set(gca, 'XTick', xticks, 'XTickLabel', compose('%.1f', xticks));

% 统一刻度数字字体：Times New Roman
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);

% 轴标签（t 斜体、s 正体；Y 轴英文 No）
xlabel('\it{t}\rm{, s}', 'FontName', 'Times New Roman', 'FontSize', 14, 'Interpreter','tex');
ylabel('No, -',         'FontName', 'Times New Roman', 'FontSize', 14);

% 边框加粗
box on;
set(gca, 'LineWidth', 1.8);

%% 配色与色标（解决“位移”中文乱码）
try, colormap(jet); catch, colormap(jet); end
cb = colorbar;

% 色标刻度数字：Times New Roman
set(cb, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.2);
cb.Label.String      = '\it{u}\rm{, m}';    
cb.Label.FontSize    = 14;
cb.Label.FontName    = 'Times New Roman'; 

%% 色标以 0 为中心（若有有效数据）
valid = Ce(~isnan(Ce));
if ~isempty(valid)
    v = max(abs(valid));
    if v > 0, caxis([-v, v]); end
end

