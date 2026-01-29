clear; clc; close all
% 文件名（请替换为你的实际文件名）
filename = 'point_ke.txt';

fid = fopen(filename, 'r');
fgetl(fid);  % 跳过标题行
data = textscan(fid, '%f %f');  % time 和 ke
fclose(fid);

time = data{1};
ke = data{2};

% === 采样信息 ===
dt_array = diff(time);
if max(abs(dt_array - mean(dt_array))) > 1e-6
    error('时间间隔不均匀，不能直接用FFT');
end
dt = dt_array(1);
Fs = 1 / dt;
N = length(ke);
ke_detrended = ke - mean(ke);

% === 方法1: 手动 FFT 计算 PSD ===
Y = fft(ke_detrended);
f_fft = (0:N-1)*(Fs/N);
half_N = floor(N/2)+1;
f_fft = f_fft(1:half_N);
PSD_fft = (abs(Y(1:half_N)).^2)/(N*Fs);
PSD_fft(2:end-1) = 2 * PSD_fft(2:end-1);  % 单边谱

% === 方法2: pwelch 计算 PSD ===
[PSD_welch, f_welch] = pwelch(ke_detrended, [], [], [], Fs);

% === Kolmogorov -5/3 斜率参考线 ===
% 选择合适的惯性子区频率范围
f_ref = [10, 150];  % 可根据你的数据调整
% 找对应的 PSD 起始参考值
idx_start = find(f_fft >= f_ref(1), 1, 'first');
C = PSD_fft(idx_start);  % 起点对齐
ref_line = C * (f_ref / f_ref(1)).^(-5/3);

% === 绘图 ===
figure;
loglog(f_fft, PSD_fft, 'b-', 'LineWidth', 1.2); hold on;
% loglog(f_welch, PSD_welch, 'r--', 'LineWidth', 1.2);
loglog(f_ref, ref_line, 'k-.', 'LineWidth', 1.5);  % -5/3参考线

% 标签和图例
text(f_ref(2), ref_line(2)*1.5, '\propto f^{-5/3}', ...
    'HorizontalAlignment', 'right', 'FontSize', 12);

xlabel('Frequency [Hz]');
ylabel('Power Spectral Density');
title('PSD Comparison with Kolmogorov Slope');
legend('Manual FFT',  'Kolmogorov -5/3');
grid on;
xlim([1, Fs/2]);  % 可视频率范围