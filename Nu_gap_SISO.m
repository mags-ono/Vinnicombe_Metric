%% ========================================================================
% Helper functions
%  ========================================================================

% Change format of the plant
function Pcell = tf_to_cell_struct(P)
    [b,a] = tfdata(P,'v');
    Pcell = cell(1,1);
    Pcell{1,1} = struct('b', b, 'a', a);
end

function y = apply_P(P, u, sigma)
    y = simulate_plant_general_discrete(P, u, sigma);
end

%Plant simulator
function y_time = simulate_plant_general_discrete(G0, u_time, sigma)
    [n_outputs, n_inputs] = size(G0);
    n_time_steps = size(u_time, 2);
    y_time = zeros(n_outputs, n_time_steps);

    for i = 1:n_outputs
        for j = 1:n_inputs
            b = G0{i, j}.b;
            a = G0{i, j}.a;
            resp = filter(b, a, u_time(j, :));
            y_time(i, :) = y_time(i, :) + resp;
        end
        % add measurement noise per output channel
        %y_time(i,:) = y_time(i,:) + sigma * randn(1, n_time_steps); %std
        y_time(i,:) = y_time(i,:) + sqrt(sigma) * randn(1, n_time_steps); %variance
        
    end
end

% Welch approximation of the winding number (across iterations)
function acc = specacc_init(N)
    acc.N = N;
    acc.count = 0;
    acc.Suu  = zeros(N,1);
    acc.Sy1u = zeros(N,1);
    acc.Sy2u = zeros(N,1);
end

function acc = specacc_update(acc, U, Y1, Y2)
    U  = U(:);  Y1 = Y1(:);  Y2 = Y2(:);
    acc.count = acc.count + 1;
    acc.Suu  = acc.Suu  + (U  .* conj(U));
    acc.Sy1u = acc.Sy1u + (Y1 .* conj(U));
    acc.Sy2u = acc.Sy2u + (Y2 .* conj(U));
end

function [G1hat, G2hat, Suu] = specacc_get_frf(acc, eps0)
    K = max(acc.count,1);
    Suu  = acc.Suu  / K;
    Sy1u = acc.Sy1u / K;
    Sy2u = acc.Sy2u / K;

    G1hat = Sy1u ./ max(Suu, eps0);
    G2hat = Sy2u ./ max(Suu, eps0);
end

function report = check_C_condition_from_acc(acc, tol_f, eps0, do_smooth, smooth_bins)
    if nargin < 4, do_smooth = false; end
    if nargin < 5, smooth_bins = 31; end

    [G1hat, G2hat, Suu] = specacc_get_frf(acc, eps0);

    % Frequency smoothing with kernel Hann (cheap, without FFT extra)
    if do_smooth
        h = hann_kernel(smooth_bins);
        Suu_s = conv(Suu, h, 'same');

        % smoothing :
        Num1 = conv((G1hat.*Suu), h, 'same');
        Num2 = conv((G2hat.*Suu), h, 'same');
        G1hat = Num1 ./ max(Suu_s, eps0);
        G2hat = Num2 ./ max(Suu_s, eps0);
    end

    % f1, f2
    %f1 = 1 + G1hat .* conj(G1hat);
    f1 = 1 + abs(G1hat).^2;
    f2 = 1 + G2hat .* conj(G1hat);

    f1s = fftshift(f1);
    f2s = fftshift(f2);

    report.minabs_f1 = min(abs(f1s));
    report.minabs_f2 = min(abs(f2s));
    report.wind_f1   = winding_number(f1s);
    report.wind_f2   = winding_number(f2s);

    report.inC = (report.minabs_f1 > tol_f) && (report.minabs_f2 > tol_f) ...
                 && (report.wind_f1 == report.wind_f2);
end

function w = winding_number(z)
    z = z(:);
    z = [z; z(1)];                                  % we close the contour
    dphi = unwrap(angle(z(2:end)./z(1:end-1)));     % robust increments
    w = round(sum(dphi)/(2*pi));
end

%hann smoother
function h = hann_kernel(L)
    if mod(L,2)==0, L=L+1; end
    h = hann(L);
    h = h / sum(h);
end

% Evalute discrete TF expressed using 'z' for last example
function G = frf_dt_tfdata_z(P, w)
    [b,a] = tfdata(P,'v');
    z = exp(1j*w);                 % z = e^{jw}
    num = polyval(b, z);
    den = polyval(a, z);
    G = (num ./ den).';
end

%% nu_gap_time_with_C_welch.m  (SISO, data-driven, time-domain experiments)
clear all;
clc;
rng(1);


%% Case Study: Poles in RHP
% -------------------------
% Ts = 1;
% z  = tf('z',Ts);
% P1 = tf(1,1,Ts);
% P2 = 2/z;                          % 2*z^{-1}
% 
% Ts = 0.001; z = tf('z',Ts);
% 
% P1 = (1 - 1/z) / (1 - 0.8/z);     % cero en z=1.2 (fuera), polo estable
% P2 = 1.8*(1 - 1/z) / (1 - 0.8/z) * 1/z;

%% Case Study: Heavy-Duty Gas Turbine CASE A.2. (Unknown Models) General Electric Models 6F vs 9F
% 
% Ts=0.05;
% s = tf('s');
% Kt = 1.3;            % Rowen stiff-system gain (incremental): DeltaP = 1.3*DeltaWF
% 
% % Parameters (as you requested)
% Tf56  = 0.1;   Kf56  = 1;   Tcd56 = 0.1;  ecr56=0.01;   % Model 5&6 
% Tf79  = 0.1;   Kf79  = 0;   Tcd79 = 0.2;  ecr79=0.01; % Model 7&9
% 
% % Model 5&6 (Kf56 = 1)
% Gf_open_56 = 1/(Tf56*s + 1);
% Gf_cl_56   = feedback(Gf_open_56, Kf56, -1);   % G/(1 + Kf*G)  (negative feedback)
% Gcd_56 = 1/(Tcd56*s + 1);
% 
% Pc56 = Kt * Gf_cl_56 * Gcd_56;
% 
% % Model 7&9 (Kf = 0)
% Gf_79  = 1/(Tf79*s + 1);                   % since Kf79=0, no feedback needed
% Gcd_79 = 1/(Tcd79*s + 1);
% 
% Pc79 = Kt * Gf_79 * Gcd_79;                % u -> DeltaP
% 
% % Discretize (Ts = 1s)
% 
% % P1 = c2d(Pc56, Ts);
% % P2 = c2d(Pc79, Ts);
% P2 = c2d(Pc56, Ts);
% P1 = c2d(Pc79, Ts);

%% Case Study: Heavy-Duty Gas Turbine CASE A.1. (Nominal Model) General Electric Tavakoli Parameters vs Khormali Parameters
Ts = 0.05;
s = tf('s');

% Parameters (as you requested)

Tf  = 0.26;     Kf  = 0;   Tcd = 0.16;      Kt=1.158;    ecr=0.005;   % Real Plant
Tf0 = 0.41;     Kf0  = 0;  Tcd0 = 0.0823;  Kt0=0.8302;   ecr0=0.005;  % Nominal Plant

% Real Model
Gf = 1/(Tf*s + 1);
Gcd = 1/(Tcd*s + 1);
Pc = Kt * Gf * Gcd;

% Nominal Model
Gf0  = 1/(Tf0*s + 1);                  
Gcd0 = 1/(Tcd0*s + 1);

Pc0 = Kt0 * Gf0 * Gcd0;               

% Discretize (Ts = 1s)

P1 = c2d(Pc, Ts, 'zoh');
P2 = c2d(Pc0, Ts, 'zoh');

%%
% True nu-gap (toolbox baseline)
[~, nugap_true] = gapmetric(P1,P2);
fprintf('gapmetric(P1,P2) = %.6f\n', nugap_true);

% Convert tf -> cell{1,1} struct with fields b,a (for filter)
P1c = tf_to_cell_struct(P1);
P2c = tf_to_cell_struct(P2);

% -------------------------
% Parameters
M        = 150;     % #iterations/experiments
N        = 9000;    % time-signal length --- 1000 for last, 9000 the rest
gamma    = 1;       % input-energy scaling
sigma_y1 = 0.01;    % noise std on plant 1 output  (<<< noise ON)
sigma_y2 = 0.01;    % noise std on plant 2 output  (<<< noise ON)

% Spectro
wfull = 2*pi*(0:N-1)'/N;          % grilla correcta del FFT completo: [0,2pi)
kpos  = 1:(floor(N/2)+1);         % bins de frecuencia no negativa: [0,pi]

eps0 = 1e-12;                          % avoid 0 division
L = sqrt(gamma) * sqrt(N);             % target ||u||2

% -------------------------
% Init u0 (time domain, real)
u = randn(1,N);
u = (L / max(norm(u,2),eps0)) * u;

delta_hat = zeros(M,1);

%
inputs = zeros(N,M);



% ============================================================
% Welch-manual (across iterations) + C-condition check
% ============================================================
acc  = specacc_init(N);

Nacc = 10;
tol_f = 1e-2;             % threshold for min|f|
smooth_bins = 31;         % kernel Hann in bins (15/31/51)

% -------------------------
% Monte Carlo parameters
Nmc = 100;   % number of Monte Carlo runs

% All trajectories saved
delta_hat_mc = zeros(M, Nmc);

% We save the first trajectory for peak analysis
delta_hat_single = zeros(M,1);
inputs_single    = zeros(N,M);
U_single         = [];
Y2_single        = [];

% ============================================================
% Monte Carlo loop 
% ============================================================
for mc = 1:Nmc
    
    rng(mc);   % distinta seed en cada corrida
    
    % -------------------------
    % Init u0 (time domain, real)
    u = randn(1,N);
    u = (L / max(norm(u,2),eps0)) * u;

    delta_hat_run = zeros(M,1);
    inputs_run    = zeros(N,M);

    % Welch-manual + C-condition check (reiniciar en cada corrida)
    acc  = specacc_init(N);

    % ============================================================
    % Main loop
    % ============================================================
    for n = 1:M
        % Experiments on time
        y1 = apply_P(P1c, u, sigma_y1);
        y2 = apply_P(P2c, u, sigma_y2);
        y1 = y1(1,:);  % SISO
        y2 = y2(1,:);

        % DFT
        U  = fft(u);
        Y1 = fft(y1);
        Y2 = fft(y2);

        % -------- Welch manual: acumular espectros (sin FFT extra)
        if n <= Nacc
            acc = specacc_update(acc, U, Y1, Y2);
        end

        % -------- C-condition check
        if n==Nacc
            rep = check_C_condition_from_acc(acc, tol_f, eps0, true, smooth_bins);
            fprintf('MC=%d | C-check n=%d: min|f1|=%.2e min|f2|=%.2e | wind1=%d wind2=%d | inC=%d\n', ...
                mc, n, rep.minabs_f1, rep.minabs_f2, rep.wind_f1, rep.wind_f2, rep.inC);

            if rep.inC==0
                delta_hat_run(n:M,1) = 1;
                break
            end
        end

        % -------- nu-gap update
        denom   = sqrt(abs(U).^2 + abs(Y1).^2 + eps0) .* sqrt(abs(U).^2 + abs(Y2).^2 + eps0);
        U_tilde = (abs(U).^2) .* (Y1 - Y2) ./ denom;

        % For input distribution
        U_n = U_tilde * N / max(norm(U_tilde,2),eps0);
        inputs_run(:,n) = abs(U_n).^2;

        % IDFT
        u_tilde = ifft(U_tilde, 'symmetric');

        % store estimate
        delta_hat_run(n) = norm(u_tilde,2) / L;

        % normalize for next iteration
        u = (L / max(norm(u_tilde,2),eps0)) * u_tilde(:).';
    end

    % Save trajectory
    delta_hat_mc(:,mc) = delta_hat_run;

    % Save first iteration
    if mc == 1
        delta_hat_single = delta_hat_run;
        inputs_single    = inputs_run;
        U_single         = U;
        Y2_single        = Y2;
    end
end

% MC AVERAGE
delta_hat_mean = mean(delta_hat_mc, 2);

% -------------------------
% Plot convergence: single run + Monte Carlo mean + MATLAB
fig = figure;
set(fig, 'Units', 'inches', 'Position', [1 1 3.0 2.2]);
hold on; box on;

h1 = plot(delta_hat_single, 'k--', 'LineWidth', 0.5, ...
    'DisplayName', '$\hat{\delta_\nu}$ Estimate ');

h2 = plot(delta_hat_mean, 'b', 'LineWidth', 0.6, ...
    'DisplayName', sprintf('Mean of %d MC', Nmc));

h3 = yline(nugap_true, 'r--', 'LineWidth', 0.8, ...
    'DisplayName', '$\delta_\nu$ MATLAB');

% Ejes
ax = gca;
ax.FontSize = 4;
ax.LineWidth = 0.1;
ax.Position = [0.22 0.22 0.60 0.64];

xlabel('Iteration $n$', 'Interpreter', 'latex', 'FontSize', 5);
ylabel('$\hat{\delta}_\nu$', 'Interpreter', 'latex', 'FontSize', 5);

ax.XRuler.Label.Units = 'normalized';
ax.YRuler.Label.Units = 'normalized';
ax.XRuler.Label.Position(2) = -0.055;
ax.YRuler.Label.Position(1) = -0.085;
ax.YAxis.TickLabelGapOffset = 1;

lgd = legend([h1 h2 h3], 'Location', 'southeast', ...
    'Interpreter', 'latex', 'Box', 'on');
lgd.FontSize = 4;
lgd.LineWidth = 0.2;
lgd.ItemTokenSize = [12 10];

%Nominal
xpos = 0.4*M;

%Unknown Models
%xpos = 0.6*M;

%RHP zero
%xpos=11;

text(xpos, nugap_true + 0.0015, ...
    sprintf('$\\delta_\\nu = %.4f$', nugap_true), ...
    'Interpreter', 'latex', ...
    'FontSize', 5, ...
    'Color', 'r', ...
    'BackgroundColor', 'w', ...
    'Margin', 1, ...
    'VerticalAlignment', 'bottom');

%Nominal
ylim([0.16, 0.24]);   % ajusta según el caso
xlim([0, 110]);
xticks([0 20 40 60 80 100]); %nominal

%Poles RHP
% ylim([min(delta_hat_mean), 1.01]);
% xlim([0, 15]);
% xticks([0 2 4 6 8 10 12 14]);

%Unknown
% ylim([0.28, 0.35]); 
%xlim([0, 150]);
%xticks([0 20 40 60 80 100 120 140]);
 

% =========================
% Input distributions
% =========================
figure;
plot(wfull, inputs_single(:,10), 'LineWidth', 1.0); hold on;
plot(wfull, inputs_single(:,50), 'LineWidth', 1.0);
plot(wfull, inputs_single(:,M),  'LineWidth', 1.0);
set(gca,'FontSize',16);
title('Input distribution');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');
legend('Iteration 10', 'Iteration 50', 'Iteration 150');
xlim([0 2*pi]);

% =========================
% Peak analysis
% =========================

[peak_mag, loc_peak_pos] = max(inputs_single(kpos,M));
fp     = kpos(loc_peak_pos);     % índice del peak físico (positivo)
w_peak = wfull(fp);

if fp == 1 || (mod(N,2)==0 && fp == N/2 + 1)
    fp_mirror = fp;   % DC o Nyquist
else
    fp_mirror = N - fp + 2;
end
w_peak_mirror = wfull(fp_mirror);

if abs(U(fp)) < 1e-12
    warning('U(w_peak) too small; Y2/U is numerically unstable.');
end

p0 = Y2(fp)/U(fp);    % Y2 asumido nominal

fprintf('\n===== PEAK ANALYSIS =====\n');
fprintf('Peak input-distribution magnitude = %.6e\n', peak_mag);
fprintf('Positive peak freq  = %.6e rad/sample\n', w_peak);
fprintf('Mirror peak freq    = %.6e rad/sample\n', w_peak_mirror);
fprintf('|p0(w_peak)|        = %.6f\n', abs(p0));
fprintf('angle p0(w_peak)    = %.6f rad\n', angle(p0));

% ------ GROUND TRUTH WINDING NUMBER-----

Nw = 65536;
w  = linspace(-pi, pi, Nw+1); w(end)=[];

G1 = frf_dt_tfdata_z(P1, w);
G2 = frf_dt_tfdata_z(P2, w);

f1 = 1 + G1.*conj(G1);             % =2
f2 = 1 + G2.*conj(G1);             % = 1 + 2e^{-jw}
% 

wind1 = winding_number(f1);
wind2 = winding_number(f2);

fprintf('=====GROUND TRUTH Winding Number===== \n');
fprintf('tfdata f1: wind=%d \n', wind1);
fprintf('tfdata f2: wind=%d \n', wind2);

