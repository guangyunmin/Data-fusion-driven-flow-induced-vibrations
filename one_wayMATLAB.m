% ------------------------------------------------------------------------
clear; clc;close all
zoneIndex = 1;                % 选择第 12 个 zone
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

%% ---------- 基本参数 ----------
L_TOTAL   = 0.99;             % 梁长 [m]
BIN_COUNT = 250;              
BIN_SIZE  = L_TOTAL / BIN_COUNT;
N_ELE     = BIN_COUNT;
N_NODE    = N_ELE + 1;
NDPN      = 4;                % 每节点自由度 (v, w, θw, θv)
N_DOFG    = N_NODE * NDPN;
N_FREE    = N_DOFG - 8;       % 两端各 4 个 DOF 固定
MID_NODE  = floor((N_NODE-1)/2);

%% ---------- 材料 & 流体参数 ----------
E_MOD     = 107e9;            % 弹性模量 [Pa]
RHO_STL   = 8400;             % 钢材密度 [kg/m^3]
D_INNER   = 0;                % 管内径 [m]
D_OUTER   = 0.009;           % 管外径 [m]
Iy        = pi*(D_OUTER^4 - D_INNER^4)/64;
Iz        = Iy;
A         = pi*(D_OUTER^2 - D_INNER^2)/4;

RHO_FLUID = 998.2;             % 流体密度 [kg/m^3]
aa=(1.07+0.56*(16.7/9)*(16.7/9));
A_DISP    = pi*(D_OUTER/2)^2*(aa^2+1)/(aa^2-1);

%% ---------- Newmark 参数 ----------
GAMMA_NM = 0.5;
BETA_NM  = 0.25;

%% ---------- 自由度映射 ----------
dof2free = -ones(N_DOFG,1);   % -1 表示约束
free_dof = zeros(N_FREE,1);
k = 1;
for n = 1:N_NODE
    for loc = 0:NDPN-1
        glob = (n-1)*NDPN + loc + 1;
        if n==1 || n==N_NODE
            dof2free(glob) = -1;
        else
            dof2free(glob) = k-1;
            free_dof(k)    = glob;
            k = k+1;
        end
    end
end

%% ---------- 装配全局刚度矩阵 & 质量矩阵 ----------
Kf        = zeros(N_FREE, N_FREE);
Mf_struct = zeros(N_FREE, N_FREE);
Mf_add    = zeros(N_FREE, N_FREE);

for e = 1:N_ELE
    Le = BIN_SIZE;
    % 局部刚度矩阵 kv, kw
    kv = [12     6*Le  -12    6*Le;
          6*Le   4*Le^2 -6*Le 2*Le^2;
         -12    -6*Le   12    -6*Le;
          6*Le   2*Le^2 -6*Le 4*Le^2];
    kw = kv;
    vidx = [1 4 5 8];
    widx = [2 3 6 7];
    Ke   = zeros(8,8);
    for i = 1:4
        for j = 1:4
            Ke(vidx(i),vidx(j)) = kv(i,j)*E_MOD*Iz/Le^3;
            Ke(widx(i),widx(j)) = kw(i,j)*E_MOD*Iy/Le^3;
        end
    end

    % 结构质量矩阵
    Mbase = [156    22*Le   54    -13*Le;
             22*Le  4*Le^2 13*Le -3*Le^2;
             54     13*Le   156   -22*Le;
            -13*Le -3*Le^2 -22*Le 4*Le^2];
    Me = zeros(8,8);
    for i = 1:4
        for j = 1:4
            Me(vidx(i),vidx(j)) = Mbase(i,j)*RHO_STL*A*Le/420;
            Me(widx(i),widx(j)) = Mbase(i,j)*RHO_STL*A*Le/420;
        end
    end

    % 流体附加质量
    Me_add_loc = zeros(8,8);
    for i = 1:4
        for j = 1:4
            Me_add_loc(vidx(i),vidx(j)) = Mbase(i,j)*RHO_FLUID*A_DISP*Le/420;
            Me_add_loc(widx(i),widx(j)) = Mbase(i,j)*RHO_FLUID*A_DISP*Le/420;
        end
    end

    % 装配到全局矩阵
    idx8 = ((e-1)*4+1):((e-1)*4+8);
    for i = 1:8
        fi = dof2free(idx8(i));
        if fi<0, continue; end
        for j = 1:8
            fj = dof2free(idx8(j));
            if fj<0, continue; end
            Kf(fi+1,fj+1)        = Kf(fi+1,fj+1)        + Ke(i,j);
            Mf_struct(fi+1,fj+1) = Mf_struct(fi+1,fj+1) + Me(i,j);
            Mf_add   (fi+1,fj+1) = Mf_add   (fi+1,fj+1) + Me_add_loc(i,j);
        end
    end
end

% 总质量矩阵
Mf = Mf_struct + Mf_add;

%% ---------- Rayleigh 阻尼 ----------
zeta = 0.013;
[V,D] = eig(Kf, Mf);
w1    = sqrt(D(1,1));
w3    = sqrt(D(min(3,end),min(3,end)));
a0    = 2*zeta*w1*w3/(w1+w3);
a1    = 2*zeta/(w1+w3);
Cf    = a0*Mf + a1*Kf;

%% ---------- 读入流体力文件 ----------

t        = time;
Fy_bin   = fx(:,:,zoneIndex);
Fz_bin   = fy(:,:,zoneIndex);
nstep    = length(t);

%% ---------- 装配全局力 Q_hist ----------
Q_hist = zeros(N_FREE, nstep);
for k = 1:nstep
    Q = zeros(N_FREE,1);
    for n = 1:N_ELE
        fy = Fy_bin(k,n) / BIN_SIZE;
        fz = Fz_bin(k,n) / BIN_SIZE;

        fe_w = fz * Le/2 * [1; Le/6; 1; -Le/6];
        fe_v = fy * Le/2 * [1; Le/6; 1; -Le/6];
        Qe   = zeros(8,1);
        for i = 1:4
            Qe(widx(i)) = fe_w(i);
            Qe(vidx(i)) = fe_v(i);
        end

        for j = 1:8
            gid = dof2free((n-1)*4 + j);
            if gid>=0
                Q(gid+1) = Q(gid+1) + Qe(j);
            end
        end
    end
    Q_hist(:,k) = Q;
end

%% ---------- 静态刚度解算（第1步载荷） ----------
Q1 = Q_hist(:,1);
uu_static = Kf \ Q1;
idv_mid_static = dof2free(MID_NODE*4+1);
idw_mid_static = dof2free(MID_NODE*4+2);
disp_v_mid_stat = uu_static(idv_mid_static+1);
disp_w_mid_stat = uu_static(idw_mid_static+1);
fprintf('【静态】中点 v=%.3e m, w=%.3e m\n', disp_v_mid_stat, disp_w_mid_stat);

%% ---------- Newmark-Beta MDOF 积分 ----------
dt = t(2)-t(1);
[uu, vv, aa] = NewmarkBeta1_MDOF(Mf, Cf, Kf, Q_hist, dt, nstep);

% 提取中点响应
idv_mid = dof2free(MID_NODE*4 + 1);
idw_mid = dof2free(MID_NODE*4 + 2);
assert(idv_mid>=0, '中点 v DOF 非自由！');
assert(idw_mid>=0, '中点 w DOF 非自由！');
v_mid = uu(idv_mid+1, :);
w_mid = uu(idw_mid+1, :);

%% ---------- 时域响应绘图 ----------
f1 = w1/(2*pi);
figure;
subplot(2,1,1);
plot(t, w_mid); xlabel('Time [s]'); ylabel('w [m]');
title('中点 w 位移响应');
xl = xlim; yl = ylim;
text(xl(1)+0.05*(xl(2)-xl(1)), yl(1)+0.9*(yl(2)-yl(1)), ...
     sprintf('f_1=%.2f Hz',f1), 'FontSize',12,'FontWeight','bold');

subplot(2,1,2);
plot(t, v_mid); xlabel('Time [s]'); ylabel('v [m]');
title('中点 v 位移响应');

%% ---------- PSD 分析与绘图（FFT 手动计算） ----------
% 采样频率
fs = 1/dt;

% 信号长度
N = length(w_mid);

% 对 w_mid 做 FFT
Yw = fft(w_mid);
% 双边功率谱
P2_w = (abs(Yw)/N).^2;
% 单边功率谱
P1_w = P2_w(1:floor(N/2)+1);
P1_w(2:end-1) = 2*P1_w(2:end-1);
% 频率向量
f_w = fs*(0:floor(N/2))/N;

% 对 v_mid 做同样处理
Yv = fft(v_mid);
P2_v = (abs(Yv)/N).^2;
P1_v = P2_v(1:floor(N/2)+1);
P1_v(2:end-1) = 2*P1_v(2:end-1);
f_v = f_w;  % 采样率和长度相同，频率向量一致

% 绘图
figure;
subplot(2,1,1);
loglog(f_w, P1_w);
xlabel('Frequency [Hz]');
ylabel('PSD [dB/Hz]');
title('Midpoint w 位移 PSD (FFT 计算)');
grid on;

subplot(2,1,2);
loglog(f_v, P1_v);
xlabel('Frequency [Hz]');
ylabel('PSD [dB/Hz]');
title('Midpoint v 位移 PSD (FFT 计算)');
grid on;






% --- 基本参数 ---
Nseg = 250;
BIN_SIZE = L_TOTAL / Nseg;

% Fy_bin: Nseg x Nt，行是段，列是时间
% 你原来转置过，但 surf 不需要转置，只要 meshgrid 对齐即可
x = (BIN_SIZE/2) : BIN_SIZE : (L_TOTAL - BIN_SIZE/2);   % 1 x 250


% 确保数据尺寸正确
F = Fy_bin';   % [250 x Nt]

% 创建网格
[Tg, Xg] = meshgrid(t, x);

% --- surf 绘图 ---
figure('Position',[100 100 900 520]);
surf(Tg, Xg, F, 'EdgeColor','none');  % 去掉网格线
view(2);                              % 从上往下看，变成云图
set(gca,'YDir','normal');
axis tight;
xlabel('时间 t [s]','FontSize',12);
ylabel('梁长位置 x [m]','FontSize',12);
title('流体力密度分布','FontSize',14,'FontWeight','bold');
cb = colorbar; ylabel(cb,'力密度 [N/m]');
colormap(jet);                        % 你原本用 jet，也可以换 parula/turbo
set(gca,'FontSize',12,'LineWidth',1);

% figure
% plot(t,Fy_bin(:,150))

idx = t<=0.5;
rms_w = sqrt(mean(w_mid(idx).^2));
rms_v = sqrt(mean(v_mid(idx).^2));
fprintf('0–0.4s 内 RMS: w=%.6e m, v=%.6e m\n', rms_w, rms_v);

disp('仿真及分析完成！');



%% ---------- 1/4 跨处位移提取（Hermite 插值，修正版） ----------
xq = L_TOTAL/4;             % 1/4 跨处物理位置
Le = BIN_SIZE;              % 单元长度
e_q = min(max(floor(xq/Le)+1, 1), N_ELE);   % 含有 xq 的单元编号(1-based)
n1  = e_q;                                   % 左节点
n2  = e_q + 1;                               % 右节点
xi  = (xq - (n1-1)*Le)/Le;                   % 局部坐标 xi∈[0,1]

% Euler-Bernoulli 梁的 4 自由度 Hermite 形函数（位移场）
H1 = 1 - 3*xi^2 + 2*xi^3;
H2 = Le*(xi - 2*xi^2 + xi^3);
H3 = 3*xi^2 - 2*xi^3;
H4 = Le*(-xi^2 + xi^3);

% 取相邻两节点的自由度时间历程（端部被约束则按 0 处理）
v1   = local_get_dof_series(n1, 1, dof2free, uu, nstep, NDPN);  % v(n1)
thv1 = local_get_dof_series(n1, 4, dof2free, uu, nstep, NDPN);  % θv(n1)
v2   = local_get_dof_series(n2, 1, dof2free, uu, nstep, NDPN);  % v(n2)
thv2 = local_get_dof_series(n2, 4, dof2free, uu, nstep, NDPN);  % θv(n2)

w1_  = local_get_dof_series(n1, 2, dof2free, uu, nstep, NDPN);  % w(n1)
thw1 = local_get_dof_series(n1, 3, dof2free, uu, nstep, NDPN);  % θw(n1)
w2_  = local_get_dof_series(n2, 2, dof2free, uu, nstep, NDPN);  % w(n2)
thw2 = local_get_dof_series(n2, 3, dof2free, uu, nstep, NDPN);  % θw(n2)

% Hermite 插值得到 1/4 跨处的时程
v_quarter = H1.*v1 + H2.*thv1 + H3.*v2 + H4.*thv2;
w_quarter = H1.*w1_ + H2.*thw1 + H3.*w2_ + H4.*thw2;

% ---------- 1/4 跨处时域响应 ----------
figure;
subplot(2,1,1);
plot(t, w_quarter, 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('w @ L/4 [m]');
title('1/4 跨处 w 位移响应'); grid on;

subplot(2,1,2);
plot(t, v_quarter, 'LineWidth', 1.2);
xlabel('Time [s]'); ylabel('v @ L/4 [m]');
title('1/4 跨处 v 位移响应'); grid on;

%% ---------- 1/4 跨处 PSD（沿用你当前 FFT/单边谱/对数坐标的做法） ----------
fs_q = 1/dt;

% w @ L/4
Nq   = numel(w_quarter);
Ywq  = fft(w_quarter);
P2_wq = (abs(Ywq)/Nq).^2;                   % 双边功率谱
P1_wq = P2_wq(1:floor(Nq/2)+1);             % 单边
if Nq >= 4, P1_wq(2:end-1) = 2*P1_wq(2:end-1); end
f_wq = fs_q*(0:floor(Nq/2))/Nq;

% v @ L/4
Yvq  = fft(v_quarter);
P2_vq = (abs(Yvq)/Nq).^2;
P1_vq = P2_vq(1:floor(Nq/2)+1);
if Nq >= 4, P1_vq(2:end-1) = 2*P1_vq(2:end-1); end
f_vq = f_wq;

% 绘图：保持你“中点 PSD”同样风格（loglog + 标签）
figure;
subplot(2,1,1);
loglog(f_wq, P1_wq, 'LineWidth', 1.2);
xlabel('Frequency [Hz]');
ylabel('PSD [dB/Hz]');      % 注：这是原代码的标签习惯；数值仍是线性 PSD
title('1/4 跨处 w 位移 PSD (FFT 计算)'); grid on;

subplot(2,1,2);
loglog(f_vq, P1_vq, 'LineWidth', 1.2);
xlabel('Frequency [Hz]');
ylabel('PSD [dB/Hz]');
title('1/4 跨处 v 位移 PSD (FFT 计算)'); grid on;

% ===== 辅助函数：安全获取指定节点/分量的时间历程（若为约束DOF则返回全零） =====
function s = local_get_dof_series(node, comp, dof2free, uu, nstep, NDPN)
    idx = dof2free((node-1)*NDPN + comp);
    if idx >= 0
        s = uu(idx+1, :);
    else
        s = zeros(1, nstep);
    end
end


%% ===== Newmark-Beta 多自由度主程序 =====
function [uu, vv, aa] = NewmarkBeta1_MDOF(M, C, K, Q_hist, dt, nstep)
    beta  = 0.25; gamma = 0.5;
    NDOF  = size(M,1);
    uu = zeros(NDOF,nstep);
    vv = zeros(NDOF,nstep);
    aa = zeros(NDOF,nstep);
    a0 = 1/(beta*dt^2);
    a1 = gamma/(beta*dt);
    a2 = 1/(beta*dt);
    a3 = 1/(2*beta)-1;
    a4 = gamma/beta-1;
    a5 = dt*(gamma/(2*beta)-1);
    a6 = dt*(1-gamma);
    a7 = dt*gamma;
    Keff = K + a0*M + a1*C;
    for k = 1:nstep-1
        Qeff = Q_hist(:,k+1) + ...
               M*(a0*uu(:,k)+a2*vv(:,k)+a3*aa(:,k)) + ...
               C*(a1*uu(:,k)+a4*vv(:,k)+a5*aa(:,k));
        uu(:,k+1) = Keff \ Qeff;
        aa(:,k+1) = a0*(uu(:,k+1)-uu(:,k)) - a2*vv(:,k) - a3*aa(:,k);
        vv(:,k+1) = vv(:,k) + a6*aa(:,k) + a7*aa(:,k+1);
    end
end
