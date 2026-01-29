% ------------------------------------------------------------------------
% read_fluid_forces.m
% 读取 fluid_forces.txt，将 fx, fy 按 zone 和 bin 存入数组
% ------------------------------------------------------------------------
clear; clc;close all

%% 1. 文件与参数设置
filename  = 'fluid_forces.txt';
zoneIDs   = [11,25,39,5,19,26,36,8,15,21,40,20];   % <-- 根据你的 UDF 中 WALL_ZONE_IDS 修改
numZones  = numel(zoneIDs);
numBins   = 250;            % 与 UDF 中 BIN_COUNT 保持一致

%% 2. 用 readtable 读取整个带表头的文件
T = readtable(filename, 'FileType','text', 'Delimiter','\t');

%% 3. 提取时间向量
% 假设表头第一列叫 'time'
time = T.time;  
nT   = height(T);

%% 4. 预分配 fx, fy 数组
fx = zeros(nT, numBins, numZones);
fy = zeros(nT, numBins, numZones);

%% 5. 按 zone & bin 拆列
for k = 1:numZones
    zid = zoneIDs(k);
    for b = 1:numBins
        % UDF 输出里 bin 从 0 到 numBins-1
        bin0 = b-1;
        col_fx = sprintf('z%d_fx%d', zid, bin0);
        col_fy = sprintf('z%d_fy%d', zid, bin0);
        fx(:, b, k) = T.(col_fx);
        fy(:, b, k) = T.(col_fy);
    end
end

%% 6. 绘图示例：第 12 个 zone 的总 Fx & Fy 随时间变化
zoneIndex = 12;                % 选择第 12 个 zone
zoneID    = zoneIDs(zoneIndex);

% 计算各时间步的总力
Fx_sum = sum(fx(:,:,zoneIndex), 2);
Fy_sum = sum(fy(:,:,zoneIndex), 2);

figure;

subplot(2,1,1);
plot(time, Fx_sum, 'LineWidth', 1.2);
ylabel('Total F_x (N)');
title([ 'Zone ' num2str(zoneID) ' — F_x vs Time' ]);
grid on;

subplot(2,1,2);
plot(time, Fy_sum, 'LineWidth', 1.2);
xlabel('Time (s)');
ylabel('Total F_y (N)');
title([ 'Zone ' num2str(zoneID) ' — F_y vs Time' ]);
grid on;
