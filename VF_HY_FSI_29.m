clear; clc; close all;

%% ====================== 一、配置区 ======================
% —— 数据文件（一个 zone 的拼接矩阵；每个是 250×(Ns*Nt) ）——
lowMat_fx   = "fy_zone_29_concat_LOW.mat";     % 低保真 fx 拼接矩阵
highMat_fx  = "fy_zone_29_concat_HIGH.mat";    % 高保真 fx 拼接矩阵（真值）
lowMat_fy   = "fx_zone_29_concat_LOW.mat";     % 低保真 fy 拼接矩阵
highMat_fy  = "fx_zone_29_concat_HIGH.mat";    % 高保真 fy 拼接矩阵

% —— 流速簇信息 —— 
Ul   = 2:12;                  % LOW 的流速集合（顺序需与低保真拼接一致）
Uh   = [2 4 6 9 10 12];       % HIGH 的流速集合（顺序需与高保真拼接一致）
Uval = 9;                     % 留出的高保真流速（∈ Uh）

% —— POD/RBF 超参 —— 
energy_L      = 0.99;         % LOW 能量阈值（用于确定 r）
energy_H      = 0.99;         % HIGH 能量阈值
energy_thresh = max(energy_L, energy_H);
rbf_reg       = 1e-10;        % RBF 正则
eps_scale     = [0.5 1 2 4];  % RBF eps 倍率（基于缩放后 U 的成对距离中位数）

% —— 时轴设置 —— 
dt               = 5e-4;      % 采样步长（秒）
unstable_prefix  = 51;        % 每段开头剔除帧数（LOW/HIGH 统一）
unstable_tail    = 0;         % 结尾剔除帧数

% —— 可视化/保存 —— 
save_fig      = false;
out_dir       = "fig_out";
binID         = 125;           % 指定一个 bin 做时程/PSD 对比
sum_over_bins = true;          % 是否合成“全域力”
bin_weights   = [];            % 为空则等权；也可传入 250×1 权重

% —— 梁/材料/流体参数 —— 
L_TOTAL = 0.99;               % 梁长 [m]
BIN_COUNT = 250;              % 段数（=bin 数）
BIN_SIZE  = L_TOTAL / BIN_COUNT;
N_ELE   = BIN_COUNT;
N_NODE  = N_ELE + 1;
NDPN    = 4;                  % 每节点 DOF: [v, w, theta_w, theta_v]
N_DOFG  = N_NODE * NDPN;
N_FREE  = N_DOFG - 8;         % 两端各 4 DOF 固定

E_MOD   = 107e9;              % 弹性模量 [Pa]
RHO_STL = 8400;               % 钢密度 [kg/m^3]
D_INNER = 0.0;                % 管内径 [m]
D_OUTER = 0.009;              % 管外径 [m]
Iy      = pi*(D_OUTER^4 - D_INNER^4)/64;
Iz      = Iy;
A       = pi*(D_OUTER^2 - D_INNER^2)/4;

RHO_FLUID = 998.2;            % 流体密度 [kg/m^3]
aa        = (1.07 + 0.56*(16.7/9)*(16.7/9));
A_DISP    = pi*(D_OUTER/2)^2*(aa^2+1)/(aa^2-1);  % 等效排水截面（附加质量）

% —— Newmark 与阻尼 —— 
GAMMA_NM = 0.5; BETA_NM = 0.25;
zeta     = 0.013;             % 目标阻尼比（用 w1/w3 拟合 Rayleigh）

% —— 方向映射（如你的 fx/fy 物理含义相反，在此对调）——
map_fx_to_v = true;  % true: fx→v(侧向); false: fx→w
map_fy_to_w = true;  % true: fy→w(法向);  false: fy→v

%% ====================== 二、变保真预测（fx & fy） ======================
% 返回：F_pred(250×Nt)、F_true_raw(250×Nt)、Nt、info(含 PhiH, muH, r)
[Fx_pred, Fx_true, Nt_fx, info_fx] = predict_field_for_component( ...
    lowMat_fx, highMat_fx, Ul, Uh, Uval, unstable_prefix, unstable_tail, ...
    energy_thresh, rbf_reg, eps_scale, "fx");

[Fy_pred, Fy_true, Nt_fy, info_fy] = predict_field_for_component( ...
    lowMat_fy, highMat_fy, Ul, Uh, Uval, unstable_prefix, unstable_tail, ...
    energy_thresh, rbf_reg, eps_scale, "fy");

assert(Nt_fx==Nt_fy, 'fx 与 fy 的剪裁后时间步数应一致');
Nt = Nt_fx; t = (0:Nt-1)*dt;   % 统一时间轴（稳定段起算 t=0）

%% —— 只输出：时域力对比 + 力的 PSD 对比（bin 级 & 全域）——
% Bin 级：fx
figure('Color','w','Name',sprintf('fx: 时域对比 bin=%d, U=%g',binID,Uval));
plot(t, Fx_true(binID,:),'-','LineWidth',1.4); hold on;
plot(t, Fx_pred(binID,:),'--','LineWidth',1.4); grid on;
xlabel('Time [s]'); ylabel('f_x'); legend({'真实','预测'},'Location','best');
title(sprintf('fx@bin %d — 时域：真实 vs 预测', binID));

[f_fx, P_fx_true] = compute_psd_fft(Fx_true(binID,:), dt);
[~,   P_fx_pred]  = compute_psd_fft(Fx_pred(binID,:),  dt);
mask = f_fx>0; to_dB = @(P) 10*log10(max(P, realmin));
figure('Color','w','Name',sprintf('fx: PSD 对比 bin=%d, U=%g',binID,Uval));
semilogx(f_fx(mask), to_dB(P_fx_true(mask)),'-','LineWidth',1.4); hold on;
semilogx(f_fx(mask), to_dB(P_fx_pred(mask)),'--','LineWidth',1.4); grid on;
xlabel('Frequency [Hz]'); ylabel('PSD [dB/Hz]');
legend({'真实','预测'},'Location','best');
title(sprintf('fx@bin %d — PSD：真实 vs 预测', binID));

% Bin 级：fy
figure('Color','w','Name',sprintf('fy: 时域对比 bin=%d, U=%g',binID,Uval));
plot(t, Fy_true(binID,:),'-','LineWidth',1.4); hold on;
plot(t, Fy_pred(binID,:),'--','LineWidth',1.4); grid on;
xlabel('Time [s]'); ylabel('f_y'); legend({'真实','预测'},'Location','best');
title(sprintf('fy@bin %d — 时域：真实 vs 预测', binID));

[f_fy, P_fy_true] = compute_psd_fft(Fy_true(binID,:), dt);
[~,   P_fy_pred]  = compute_psd_fft(Fy_pred(binID,:),  dt);
mask2 = f_fy>0;
figure('Color','w','Name',sprintf('fy: PSD 对比 bin=%d, U=%g',binID,Uval));
semilogx(f_fy(mask2), to_dB(P_fy_true(mask2)),'-','LineWidth',1.4); hold on;
semilogx(f_fy(mask2), to_dB(P_fy_pred(mask2)),'--','LineWidth',1.4); grid on;
xlabel('Frequency [Hz]'); ylabel('PSD [dB/Hz]');
legend({'真实','预测'},'Location','best');
title(sprintf('fy@bin %d — PSD：真实 vs 预测', binID));

% 全域：fx+fy（如启用）
if sum_over_bins
    if isempty(bin_weights), w = ones(250,1); else, w = bin_weights(:); assert(numel(w)==250,'bin_weights 需为 250×1'); end
    F_true_global = (w.' * (Fx_true + Fy_true));     % 1×Nt
    F_pred_global = (w.' * (Fx_pred + Fy_pred));     % 1×Nt

    figure('Color','w','Name',sprintf('全域力：时域 U=%g',Uval));
    plot(t, F_true_global,'-','LineWidth',1.4); hold on;
    plot(t, F_pred_global,'--','LineWidth',1.4); grid on;
    xlabel('Time [s]'); ylabel('F_{global}');
    legend({'真实','预测'},'Location','best');
    title('全域力 — 时域：真实 vs 预测');

    [fg, Pg_true] = compute_psd_fft(F_true_global, dt);
    [~,  Pg_pred] = compute_psd_fft(F_pred_global, dt);
    maskg = fg>0;
    figure('Color','w','Name',sprintf('全域力：PSD U=%g',Uval));
    semilogx(fg(maskg), to_dB(Pg_true(maskg)),'-','LineWidth',1.4); hold on;
    semilogx(fg(maskg), to_dB(Pg_pred(maskg)),'--','LineWidth',1.4); grid on;
    xlabel('Frequency [Hz]'); ylabel('PSD [dB/Hz]');
    legend({'真实','预测'},'Location','best');
    title('全域力 — PSD：真实 vs 预测');
end

%% ====================== 三、FSI：位移时域对比（真实/预测/真实(投影)） ======================
% 1) 方向映射
if map_fx_to_v
    Fv_true = Fx_true; Fw_true = Fy_true;
    Fv_pred = Fx_pred; Fw_pred = Fy_pred;
else
    Fv_true = Fy_true; Fw_true = Fx_true;
    Fv_pred = Fy_pred; Fw_pred = Fx_pred;
end

% 2) 真实力的“子空间投影”版本（分别对 fx 和 fy 用各自的 POD 子空间）
Fx_true_proj = project_to_subspace(Fx_true, info_fx.PhiH, info_fx.muH);
Fy_true_proj = project_to_subspace(Fy_true, info_fy.PhiH, info_fy.muH);
if map_fx_to_v
    Fv_true_proj = Fx_true_proj; Fw_true_proj = Fy_true_proj;
else
    Fv_true_proj = Fy_true_proj; Fw_true_proj = Fx_true_proj;
end

% 3) 自由度映射
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

% 4) 全局矩阵与阻尼
[Kf, Mf_struct, Mf_add] = assemble_KM(E_MOD, Iy, Iz, RHO_STL, A, RHO_FLUID, A_DISP, BIN_SIZE, N_ELE, NDPN, dof2free);
Mf = Mf_struct + Mf_add;
[w1, w3] = get_w1w3(Kf, Mf);
a0 = 2*zeta*w1*w3/(w1+w3);
a1 = 2*zeta/(w1+w3);
Cf = a0*Mf + a1*Kf;

% 5) 三套载荷：真实 / 预测 / 真实(投影)
Q_hist_true = build_Q_hist(Fv_true,      Fw_true,      BIN_SIZE, N_ELE, NDPN, dof2free);
Q_hist_pred = build_Q_hist(Fv_pred,      Fw_pred,      BIN_SIZE, N_ELE, NDPN, dof2free);
Q_hist_proj = build_Q_hist(Fv_true_proj, Fw_true_proj, BIN_SIZE, N_ELE, NDPN, dof2free);

% 6) 三次 Newmark
[uu_true, ~, ~] = NewmarkBeta1_MDOF(Mf, Cf, Kf, Q_hist_true, dt, Nt);
[uu_pred, ~, ~] = NewmarkBeta1_MDOF(Mf, Cf, Kf, Q_hist_pred, dt, Nt);
[uu_proj, ~, ~] = NewmarkBeta1_MDOF(Mf, Cf, Kf, Q_hist_proj, dt, Nt);

% 7) 中点位移（v/w）对比：新增“真实(投影)力驱动”
MID_NODE = floor((N_NODE-1)/2);
[idv_mid, idw_mid] = deal(dof2free(MID_NODE*NDPN + 1), dof2free(MID_NODE*NDPN + 2));
assert(idv_mid>=0 && idw_mid>=0, '中点 DOF 非自由！');

v_mid_true = uu_true(idv_mid+1,:); w_mid_true = uu_true(idw_mid+1,:);
v_mid_pred = uu_pred(idv_mid+1,:); w_mid_pred = uu_pred(idw_mid+1,:);
v_mid_proj = uu_proj(idv_mid+1,:); w_mid_proj = uu_proj(idw_mid+1,:);

figure('Color','w','Name','中点位移 — 时域（真实/预测/真实(投影)）');
subplot(2,1,1);
plot(t, w_mid_true,'-','LineWidth',1.4); hold on;
plot(t, w_mid_pred,'--','LineWidth',1.4);
plot(t, w_mid_proj,'-.','LineWidth',1.4);
grid on; xlabel('Time [s]'); ylabel('w [m]');
legend({'真实力驱动','预测力驱动','真实(投影)力驱动'},'Location','best'); title('中点 w(t)');

subplot(2,1,2);
plot(t, v_mid_true,'-','LineWidth',1.4); hold on;
plot(t, v_mid_pred,'--','LineWidth',1.4);
plot(t, v_mid_proj,'-.','LineWidth',1.4);
grid on; xlabel('Time [s]'); ylabel('v [m]');
legend({'真实力驱动','预测力驱动','真实(投影)力驱动'},'Location','best'); title('中点 v(t)');

fprintf('\n==== 完成：U=%g 的三类力对比 + 位移三曲线（含真实力投影）====\n', Uval);

%% ====== L/4 处：位移三曲线（真实/预测/真实(投影)）======
xq = L_TOTAL/4;                 % 物理位置
Le = BIN_SIZE;
e_q = min(max(floor(xq/Le)+1, 1), N_ELE);   % 含有 xq 的单元号(1-based)
n1  = e_q;                      % 左节点
n2  = e_q + 1;                  % 右节点
xi  = (xq - (n1-1)*Le)/Le;      % 局部坐标 xi∈[0,1]

% 四个 Hermite 形函数（Euler-Bernoulli 梁，位移场）
H1 = 1 - 3*xi^2 + 2*xi^3;
H2 = Le*(xi - 2*xi^2 + xi^3);
H3 = 3*xi^2 - 2*xi^3;
H4 = Le*(-xi^2 + xi^3);

% —— 从三套位移结果中取出 n1/n2 的 (v, θv, w, θw) 时间历程，并 Hermite 插值到 x=L/4 ——
getv  = @(uu_,node) get_dof_series(node, 1, dof2free, uu_, size(uu_,2), NDPN);
getthv= @(uu_,node) get_dof_series(node, 4, dof2free, uu_, size(uu_,2), NDPN);
getw  = @(uu_,node) get_dof_series(node, 2, dof2free, uu_, size(uu_,2), NDPN);
getthw= @(uu_,node) get_dof_series(node, 3, dof2free, uu_, size(uu_,2), NDPN);

% 真实
v1 = getv(uu_true,n1); v2 = getv(uu_true,n2);
tv1= getthv(uu_true,n1); tv2= getthv(uu_true,n2);
w1 = getw(uu_true,n1); w2 = getw(uu_true,n2);
tw1= getthw(uu_true,n1); tw2= getthw(uu_true,n2);
v_L4_true = H1.*v1 + H2.*tv1 + H3.*v2 + H4.*tv2;
w_L4_true = H1.*w1 + H2.*tw1 + H3.*w2 + H4.*tw2;

% 预测
v1 = getv(uu_pred,n1); v2 = getv(uu_pred,n2);
tv1= getthv(uu_pred,n1); tv2= getthv(uu_pred,n2);
w1 = getw(uu_pred,n1); w2 = getw(uu_pred,n2);
tw1= getthw(uu_pred,n1); tw2= getthw(uu_pred,n2);
v_L4_pred = H1.*v1 + H2.*tv1 + H3.*v2 + H4.*tv2;
w_L4_pred = H1.*w1 + H2.*tw1 + H3.*w2 + H4.*tw2;

% 真实(投影)
v1 = getv(uu_proj,n1); v2 = getv(uu_proj,n2);
tv1= getthv(uu_proj,n1); tv2= getthv(uu_proj,n2);
w1 = getw(uu_proj,n1); w2 = getw(uu_proj,n2);
tw1= getthw(uu_proj,n1); tw2= getthw(uu_proj,n2);
v_L4_proj = H1.*v1 + H2.*tv1 + H3.*v2 + H4.*tv2;
w_L4_proj = H1.*w1 + H2.*tw1 + H3.*w2 + H4.*tw2;

% —— 绘图（只 plot，不保存）——
figure('Color','w','Name','L/4 位移 — 时域（真实/预测/真实(投影)）');
subplot(2,1,1);
plot(t, w_L4_true,'-','LineWidth',1.4); hold on;
plot(t, w_L4_pred,'--','LineWidth',1.4);
plot(t, w_L4_proj,'-.','LineWidth',1.4);
grid on; xlabel('Time [s]'); ylabel('w @ L/4 [m]');
legend({'真实力驱动','预测力驱动','真实(投影)力驱动'},'Location','best');
title('L/4 处 w(t)');

subplot(2,1,2);
plot(t, v_L4_true,'-','LineWidth',1.4); hold on;
plot(t, v_L4_pred,'--','LineWidth',1.4);
plot(t, v_L4_proj,'-.','LineWidth',1.4);
grid on; xlabel('Time [s]'); ylabel('v @ L/4 [m]');
legend({'真实力驱动','预测力驱动','真实(投影)力驱动'},'Location','best');
title('L/4 处 v(t)');

%% ====== L/4 处：力三曲线（真实/预测/真实(投影)）======
% 使用与装配一致的方向映射（侧向 v、法向 w）
b_L4 = e_q;   % L/4 所在的 bin = 对应单元号
fv_true  = Fv_true(b_L4,:);  fw_true  = Fw_true(b_L4,:);
fv_pred  = Fv_pred(b_L4,:);  fw_pred  = Fw_pred(b_L4,:);
fv_proj  = Fv_true_proj(b_L4,:);  fw_proj  = Fw_true_proj(b_L4,:);

figure('Color','w','Name','L/4 力 — 时域（真实/预测/真实(投影)）');
subplot(2,1,1);
plot(t, fv_true,'-','LineWidth',1.4); hold on;
plot(t, fv_pred,'--','LineWidth',1.4);
plot(t, fv_proj,'-.','LineWidth',1.4);
grid on; xlabel('Time [s]'); ylabel('f_v @ L/4 [N]');
legend({'真实','预测','真实(投影)'},'Location','best');
title('L/4 处 f_v(t)');

subplot(2,1,2);
plot(t, fw_true,'-','LineWidth',1.4); hold on;
plot(t, fw_pred,'--','LineWidth',1.4);
plot(t, fw_proj,'-.','LineWidth',1.4);
grid on; xlabel('Time [s]'); ylabel('f_w @ L/4 [N]');
legend({'真实','预测','真实(投影)'},'Location','best');
title('L/4 处 f_w(t)');


%% ====== L/2 & L/4 处位移 PSD（v 与 w；真实/预测/真实(投影)）======
to_dB = @(P) 10*log10(max(P, realmin));

% ---------------- L/2（中点）位移 PSD ----------------
% w @ L/2
[f_m, Pw_true_m] = compute_psd_fft(w_mid_true, dt);
[~,   Pw_pred_m] = compute_psd_fft(w_mid_pred, dt);
[~,   Pw_proj_m] = compute_psd_fft(w_mid_proj, dt);
mask_m = f_m > 0;

% v @ L/2
[~, Pv_true_m] = compute_psd_fft(v_mid_true, dt);
[~, Pv_pred_m] = compute_psd_fft(v_mid_pred, dt);
[~, Pv_proj_m] = compute_psd_fft(v_mid_proj, dt);

figure('Color','w','Name','L/2 位移 PSD（真实/预测/真实(投影)）');
subplot(2,1,1);
semilogx(f_m(mask_m), to_dB(Pw_true_m(mask_m)),'-','LineWidth',1.4); hold on;
semilogx(f_m(mask_m), to_dB(Pw_pred_m(mask_m)),'--','LineWidth',1.4);
semilogx(f_m(mask_m), to_dB(Pw_proj_m(mask_m)),'-.','LineWidth',1.4);
grid on; xlabel('Frequency [Hz]'); ylabel('PSD of w (dB/Hz)');
legend({'真实力','预测力','真实(投影)力'},'Location','best');
title('L/2 处 w 的 PSD');

subplot(2,1,2);
semilogx(f_m(mask_m), to_dB(Pv_true_m(mask_m)),'-','LineWidth',1.4); hold on;
semilogx(f_m(mask_m), to_dB(Pv_pred_m(mask_m)),'--','LineWidth',1.4);
semilogx(f_m(mask_m), to_dB(Pv_proj_m(mask_m)),'-.','LineWidth',1.4);
grid on; xlabel('Frequency [Hz]'); ylabel('PSD of v (dB/Hz)');
legend({'真实力','预测力','真实(投影)力'},'Location','best');
title('L/2 处 v 的 PSD');

% ---------------- L/4 位移 PSD ----------------
% 已有的 L/4 位移时程：w_L4_true/pred/proj、v_L4_true/pred/proj
[f_q, Pw_true_q] = compute_psd_fft(w_L4_true, dt);
[~,   Pw_pred_q] = compute_psd_fft(w_L4_pred, dt);
[~,   Pw_proj_q] = compute_psd_fft(w_L4_proj, dt);
mask_q = f_q > 0;

[~, Pv_true_q] = compute_psd_fft(v_L4_true, dt);
[~, Pv_pred_q] = compute_psd_fft(v_L4_pred, dt);
[~, Pv_proj_q] = compute_psd_fft(v_L4_proj, dt);

figure('Color','w','Name','L/4 位移 PSD（真实/预测/真实(投影)）');
subplot(2,1,1);
semilogx(f_q(mask_q), to_dB(Pw_true_q(mask_q)),'-','LineWidth',1.4); hold on;
semilogx(f_q(mask_q), to_dB(Pw_pred_q(mask_q)),'--','LineWidth',1.4);
semilogx(f_q(mask_q), to_dB(Pw_proj_q(mask_q)),'-.','LineWidth',1.4);
grid on; xlabel('Frequency [Hz]'); ylabel('PSD of w (dB/Hz)');
legend({'真实力','预测力','真实(投影)力'},'Location','best');
title('L/4 处 w 的 PSD');

subplot(2,1,2);
semilogx(f_q(mask_q), to_dB(Pv_true_q(mask_q)),'-','LineWidth',1.4); hold on;
semilogx(f_q(mask_q), to_dB(Pv_pred_q(mask_q)),'--','LineWidth',1.4);
semilogx(f_q(mask_q), to_dB(Pv_proj_q(mask_q)),'-.','LineWidth',1.4);
grid on; xlabel('Frequency [Hz]'); ylabel('PSD of v (dB/Hz)');
legend({'真实力','预测力','真实(投影)力'},'Location','best');
title('L/4 处 v 的 PSD');


%% ====================== 附：函数区 ======================
function [F_pred, F_true_raw, NtH2, info] = predict_field_for_component(lowMatFile, highMatFile, Ul, Uh, U_val, unstable_prefix, unstable_tail, energy_thresh, rbf_reg, eps_scale, tag)
    % 读取
    assert(isfile(lowMatFile),  '未找到 %s', lowMatFile);
    assert(isfile(highMatFile), '未找到 %s', highMatFile);
    ML_big = load_one_matrix(lowMatFile);    % 250×(NsL*NtL)
    MH_big = load_one_matrix(highMatFile);   % 250×(NsH*NtH)
    assert(size(ML_big,1)==250 && size(MH_big,1)==250, '行数必须为 250');

    % 基本维度
    NsL = numel(Ul); assert(mod(size(ML_big,2), NsL)==0, 'LOW 列数无法被 NsL 整除');
    NsH = numel(Uh); assert(mod(size(MH_big,2), NsH)==0, 'HIGH 列数无法被 NsH 整除');
    NtL = size(ML_big,2)/NsL; 
    NtH = size(MH_big,2)/NsH;

    % 剪裁并拉平
    ML_flat = flatten_blocks(ML_big, NtL, unstable_prefix, unstable_tail); % (250*NtL_eff)×NsL
    MH_flat = flatten_blocks(MH_big, NtH, unstable_prefix, unstable_tail); % (250*NtH_eff)×NsH
    NtL_eff = NtL - unstable_prefix - unstable_tail; 
    NtH_eff = NtH - unstable_prefix - unstable_tail; 

    % —— 从 HIGH 原始大矩阵中切出“真值 raw 段”（不投影，仅剪裁）——
    assert(ismember(U_val, Uh), 'U_val 必须属于 Uh');
    iH_val  = find(Uh==U_val,1);
    s_raw = (iH_val-1)*NtH + 1 + unstable_prefix;
    e_raw =  iH_val   *NtH - unstable_tail;
    F_true_raw = MH_big(:, s_raw:e_raw);     % 250×NtH_eff（作为“真实力”）
    NtH2 = NtH_eff;

    % —— HIGH 训练/验证划分（用于预测）——
    idxH_tr = setdiff(1:NsH, iH_val);
    Utr     = Uh(idxH_tr);         % 训练流速
    MH_tr   = MH_flat(:, idxH_tr); % 高保真训练集
    xH_val  = MH_flat(:, iH_val);  % 留出高保真（仅用于获取子空间的均值尺度）

    % —— 样本型 POD：LOW 用全部；HIGH 仅训练子集 —— 
    muL = mean(ML_flat,2); XLc = ML_flat - muL; [UL_s,SL_s,~] = svd(XLc,'econ'); sL = diag(SL_s);
    efracL = cumsum(sL.^2)/sum(sL.^2); rL = find(efracL>=energy_thresh,1,'first'); if isempty(rL), rL=numel(sL); end

    muH = mean(MH_tr,2); XHc = MH_tr - muH; [UH_s,SH_s,~] = svd(XHc,'econ'); sH = diag(SH_s);
    efracH = cumsum(sH.^2)/sum(sH.^2); rH = find(efracH>=energy_thresh,1,'first'); if isempty(rH), rH=numel(sH); end

    nL = size(UL_s,2); nH = size(UH_s,2);
    r  = min([max(rL,rH), nL, nH]);                   % 统一 r
    PhiL = UL_s(:,1:r); %#ok<NASGU>
    PhiH = UH_s(:,1:r);

    % —— 系数 —— 
    AL   = PhiL.' * (ML_flat - muL);      % r×NsL
    BH_tr= PhiH.' * (MH_tr   - muH);      % r×(NsH-1)

    % —— RBF 训练：B-A 残差（逐模态），U 归一化 —— 
    [~, idxL_tr] = ismember(Utr, Ul);
    rbfs = cell(r,1);
    U_scaler = fit_minmax(Utr);
    Utr_s = apply_minmax(U_scaler, Utr);
    eps_list = auto_eps_list(Utr_s, eps_scale);
    for m=1:r
        y_tr   = (BH_tr(m,:).' - AL(m, idxL_tr).');
        y_scaler = fit_zscore(y_tr);
        y_tr_s   = apply_zscore(y_scaler, y_tr);
        best = struct('err',inf,'eps',NaN,'alpha',[]);
        for eps = eps_list
            D = abs(Utr_s(:)-Utr_s(:)'); K = exp(-(eps*D).^2);
            alpha = (K + rbf_reg*eye(numel(Utr_s))) \ y_tr_s;
            err = norm(K*alpha - y_tr_s)/max(1e-12, norm(y_tr_s));
            if err < best.err, best = struct('err',err,'eps',eps,'alpha',alpha); end
        end
        rbfs{m} = struct('U_scaler',U_scaler,'y_scaler',y_scaler,'Utr_s',Utr_s(:),'alpha',best.alpha,'eps',best.eps);
    end

    % —— Uval 的系数预测 —— 
    [~, iL_val] = ismember(U_val, Ul); A_val = AL(:, iL_val); Rhat = zeros(r,1);
    for m=1:r, Rhat(m) = rbf_predict_norm(rbfs{m}, U_val); end
    B_pred = A_val + Rhat;

    % —— 场重构（预测）—— 
    Lh = size(PhiH,1); NtH2_chk = round(Lh/250); assert(NtH2_chk==NtH2, '时间长度不一致');
    F_pred_flat = muH + PhiH * B_pred; 
    F_pred = reshape(F_pred_flat, 250, NtH2);

    % —— 返回信息（用于后续“真实力子空间投影”）——
    info = struct('r',r,'Nt',NtH2,'tag',tag,'PhiH',PhiH,'muH',muH);
end

%% ============= FSI 装配/积分/后处理辅助 =============
function [Kf, Mf_struct, Mf_add] = assemble_KM(E_MOD, Iy, Iz, RHO_STL, A, RHO_FLUID, A_DISP, Le, N_ELE, NDPN, dof2free)
    Kf = zeros(sum(dof2free>=0)); Mf_struct = Kf; Mf_add = Kf;
    kv = [12 6*Le -12 6*Le; 6*Le 4*Le^2 -6*Le 2*Le^2; -12 -6*Le 12 -6*Le; 6*Le 2*Le^2 -6*Le 4*Le^2];
    kw = kv; vidx = [1 4 5 8]; widx = [2 3 6 7];
    for e=1:N_ELE
        Ke = zeros(8,8); Mbase = [156 22*Le 54 -13*Le; 22*Le 4*Le^2 13*Le -3*Le^2; 54 13*Le 156 -22*Le; -13*Le -3*Le^2 -22*Le 4*Le^2]; 
        Me = zeros(8,8); Me_add = zeros(8,8);
        for i=1:4
            for j=1:4
                Ke(vidx(i),vidx(j)) = kv(i,j)*E_MOD*Iz/Le^3; Ke(widx(i),widx(j)) = kw(i,j)*E_MOD*Iy/Le^3;
                Me(vidx(i),vidx(j)) = Mbase(i,j)*RHO_STL*A*Le/420; Me(widx(i),widx(j)) = Mbase(i,j)*RHO_STL*A*Le/420;
                Me_add(vidx(i),vidx(j)) = Mbase(i,j)*RHO_FLUID*A_DISP*Le/420; Me_add(widx(i),widx(j)) = Mbase(i,j)*RHO_FLUID*A_DISP*Le/420;
            end
        end
        idx8 = ((e-1)*4+1):((e-1)*4+8);
        for i=1:8
            fi = dof2free(idx8(i)); if fi<0, continue; end
            for j=1:8
                fj = dof2free(idx8(j)); if fj<0, continue; end
                Kf(fi+1,fj+1)        = Kf(fi+1,fj+1)        + Ke(i,j);
                Mf_struct(fi+1,fj+1) = Mf_struct(fi+1,fj+1) + Me(i,j);
                Mf_add(fi+1,fj+1)    = Mf_add(fi+1,fj+1)    + Me_add(i,j);
            end
        end
    end
end

function [w1, w3] = get_w1w3(Kf, Mf)
    [V,D] = eig(Kf,Mf); lam = real(diag(D)); lam = lam(lam>0); lam = sort(lam,'ascend');
    if numel(lam)<3
        w1 = sqrt(lam(1)); w3 = sqrt(lam(min(1,end)));
    else
        w1 = sqrt(lam(1)); w3 = sqrt(lam(3));
    end
end

function Q_hist = build_Q_hist(Fv_bin, Fw_bin, Le, N_ELE, NDPN, dof2free)
    % Fv_bin/Fw_bin: 250×Nt 的“每 bin 力”（单位 N/段）——将其转成分布载荷 q=N/m
    [nbin, Nt] = size(Fv_bin); assert(nbin==N_ELE, 'bin 数需等于单元数');
    Q_hist = zeros(sum(dof2free>=0), Nt);
    vidx = [1 4 5 8]; widx = [2 3 6 7];
    for k=1:Nt
        Q = zeros(sum(dof2free>=0),1);
        for e=1:N_ELE
            qv = Fv_bin(e,k) / Le;  % N/m
            qw = Fw_bin(e,k) / Le;  % N/m
            fe_v = qv * Le/2 * [1; Le/6; 1; -Le/6];
            fe_w = qw * Le/2 * [1; Le/6; 1; -Le/6];
            Qe = zeros(8,1); Qe(vidx) = fe_v; Qe(widx) = fe_w;
            idx8 = ((e-1)*4+1):((e-1)*4+8);
            for j=1:8
                gid = dof2free(idx8(j)); if gid>=0, Q(gid+1) = Q(gid+1) + Qe(j); end
            end
        end
        Q_hist(:,k) = Q;
    end
end

function [uu, vv, aa] = NewmarkBeta1_MDOF(M, C, K, Q_hist, dt, nstep)
    beta=0.25; gamma=0.5; NDOF=size(M,1);
    uu=zeros(NDOF,nstep); vv=zeros(NDOF,nstep); aa=zeros(NDOF,nstep);
    a0=1/(beta*dt^2); a1=gamma/(beta*dt); a2=1/(beta*dt); a3=1/(2*beta)-1; a4=gamma/beta-1; a5=dt*(gamma/(2*beta)-1); a6=dt*(1-gamma); a7=dt*gamma;
    Keff=K + a0*M + a1*C;
    for k=1:nstep-1
        Qeff = Q_hist(:,k+1) + M*(a0*uu(:,k)+a2*vv(:,k)+a3*aa(:,k)) + C*(a1*uu(:,k)+a4*vv(:,k)+a5*aa(:,k));
        uu(:,k+1) = Keff \ Qeff;
        aa(:,k+1) = a0*(uu(:,k+1)-uu(:,k)) - a2*vv(:,k) - a3*aa(:,k);
        vv(:,k+1) = vv(:,k) + a6*aa(:,k) + a7*aa(:,k+1);
    end
end

%% ============= PSD / I-O / 归一化工具 =============
function [f, P1] = compute_psd_fft(y, dt)
    y = y(:); N = numel(y); fs=1/dt; y = y - mean(y);
    X = fft(y); Sxx_two = (abs(X).^2) / (N*fs);
    if mod(N,2)==0
        k_max = N/2 + 1; Sxx_pos = Sxx_two(1:k_max); Sxx_pos(2:end-1) = 2*Sxx_pos(2:end-1); f = (0:(k_max-1))' * (fs/N);
    else
        k_max = (N+1)/2; Sxx_pos = Sxx_two(1:k_max); Sxx_pos(2:end)   = 2*Sxx_pos(2:end);   f = (0:(k_max-1))' * (fs/N);
    end
    P1 = Sxx_pos;
end

function M = load_one_matrix(fname)
    S = load(fname); cand={'M_zone','M','ML','MH','data'};
    for i=1:numel(cand)
        if isfield(S,cand{i}) && isnumeric(S.(cand{i})), M=S.(cand{i}); return; end
    end
    f=fieldnames(S);
    for i=1:numel(f)
        X=S.(f{i});
        if isnumeric(X) && ismatrix(X) && size(X,1)==250, M=X; return;
        end
    end
    error('文件 %s 未找到 250×N 矩阵。', fname);
end

function Xflat = flatten_blocks(Mbig, Nt, skip_head, skip_tail)
    if nargin<3, skip_head=0; end; if nargin<4, skip_tail=0; end
    [nrow, ncol] = size(Mbig); assert(nrow==250 && mod(ncol,Nt)==0, '需 250×(Ns*Nt)');
    Ns=ncol/Nt; Nt_eff = Nt - skip_head - skip_tail; assert(Nt_eff>0,'剪裁后 Nt_eff>0');
    Xflat = zeros(250*Nt_eff, Ns);
    s=1; for i=1:Ns
        e=s+Nt-1; blk=Mbig(:,s:e); idx=(1+skip_head):(Nt - skip_tail); blk=blk(:,idx);
        Xflat(:,i)=blk(:); s=e+1;
    end
end

function scaler = fit_minmax(U)
    U=U(:); scaler.min=min(U); scaler.max=max(U); if scaler.max==scaler.min, scaler.max=scaler.min+1; end
end
function Us = apply_minmax(scaler, U)
    Us = (U(:)-scaler.min) / (scaler.max - scaler.min);
end
function scaler = fit_zscore(y)
    y=y(:); scaler.mu=mean(y); scaler.sig=std(y); if scaler.sig==0 || ~isfinite(scaler.sig), scaler.sig=1; end
end
function ys = apply_zscore(scaler, y)
    ys = (y(:)-scaler.mu) / scaler.sig;
end
function y = invert_zscore(scaler, ys)
    y = ys(:)*scaler.sig + scaler.mu;
end
function eps_list = auto_eps_list(U_scaled, scales)
    U_scaled=U_scaled(:); n=numel(U_scaled);
    if n>=2
        D=abs(U_scaled - U_scaled.'); d=D(triu(true(n),1)); dmed=median(d); if dmed<=0 || ~isfinite(dmed), dmed=1; end
    else
        dmed=1;
    end
    eps0=1/dmed; eps_list = eps0*scales(:).';
end
function yq = rbf_predict_norm(rbf, Uq)
    Uq_s = apply_minmax(rbf.U_scaler, Uq);
    Dq = abs(Uq_s(:) - rbf.Utr_s(:)'); Kq = exp(-(rbf.eps * Dq).^2);
    yq_s = Kq * rbf.alpha;                    % 预测（标准化）
    yq   = invert_zscore(rbf.y_scaler, yq_s); % 反归一化
end

function F_proj = project_to_subspace(F_raw, PhiH, muH)
% 将 250×Nt 的“真实 raw 力场”正交投影到预测所用 POD 子空间 span(PhiH)
% 输入：
%   F_raw : 250×NtH_eff
%   PhiH  : (250*NtH_eff)×r  —— 高保真训练得到的基（与预测一致）
%   muH   : (250*NtH_eff)×1  —— 对应均值（与预测一致）
% 输出：
%   F_proj: 250×NtH_eff      —— 投影后的力场（与预测/真实长度一致）
    z = F_raw(:);
    coeff = PhiH' * (z - muH);
    zproj = muH + PhiH * coeff;
    Nt = size(F_raw,2);
    F_proj = reshape(zproj, 250, Nt);
end
% ===== 辅助：安全获取某节点、某分量的时间历程（若为约束 DOF 则返回全零） =====
function s = get_dof_series(node, comp, dof2free, uu, nstep, NDPN)
    % comp: 1=v, 2=w, 3=theta_w, 4=theta_v
    idx = dof2free((node-1)*NDPN + comp);
    if idx >= 0
        s = uu(idx+1, :);
    else
        s = zeros(1, nstep);
    end
end
