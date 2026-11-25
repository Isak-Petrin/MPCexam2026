%% Problem 6.2 - Static and dynamic Kalman filters, Problem 4 model, no disturbance steps
% Identified transfer function model from Problem 4.4 (+10 % steps):
%
%   G11 = K11/(T1 s+1)              G12 = K12/((T1 s+1)(T3 s+1))
%   G21 = K21/((T2 s+1)(T4 s+1))    G22 = K22/(T2 s+1)
%
% Each SISO transfer function is realized in observer canonical form
% (Problem 4.6) and the four realizations are stacked into one MIMO model
% with x = [x1; x2; x3; x4], x1, x4 in R, x2, x3 in R^2, i.e. x in R^6:
%
%   A = blkdiag(A1,A2,A3,A4),  B = [B1 0; 0 B2; B3 0; 0 B4],
%   C = [C1 C2 0 0; 0 0 C3 C4],  D = 0
%
% Disturbance input. The step experiments only identify u -> y. F3 enters
% tank 3 exactly as (1-gamma2)*F2 does, and F4 enters tank 4 as
% (1-gamma1)*F1. Hence F3 -> h1 = G12/(1-gamma2) and F4 -> h2 = G21/(1-gamma1),
% which share the states of G12 and G21:
%
%   E = [0 0; B2/(1-gamma2) 0; 0 B3/(1-gamma1); 0 0]
%
% Deviation variables: u = F12 - us with us = [300; 300] (held constant,
% du = 0), d = F34 - ds with ds = [100; 80], y = h12 - hs (+ noise).
%
% Two stochastic disturbance models (as for the Problem 5 model):
%   1) piecewise constant:  x_{k+1} = Ad x_k + Bd u_k + Ed dd_k, dd_k ~ N(0,Rdd)
%                           -> G = Ed, Q = Rdd, Qpc = Ed Rdd Ed'
%   2) SDE: dx = (A x + B u) dt + Gc dw,  Gc = E diag(sigma)
%                           -> G = I,  Q = Qsde = int_0^Ts e^{As} Gc Gc' e^{A's} ds
%   y_k = C x_k + v_k, v_k ~ N(0,Rvv)
%
% The plant is the same linear 6-state model. The states of the realization
% are not physical; only the levels h1, h2 are compared.
%
% Uses KalmanFilterDesign, KalmanFilterDynamic, KalmanFilterStatic (../)
% and DiscretizeSDE (problem_5).

clc
clear
close all

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..');
addpath(fullfile(here,'..'));                          % Kalman filter functions
addpath(fullfile(root,'problem_2','problem_2_1'));     % parameters
addpath(fullfile(root,'problem_5'));                   % DiscretizeSDE

QuadrupleTankParameters                                % gamma1, gamma2, us, ds, Rdd, Rvv, sigma, seed
load(fullfile(root,'problem_4','data','tf_identification.mat'), 'Kid', 'Tid');
load(fullfile(root,'problem_5','data','linear_models.mat'), 'Ts', 'hs');
P5 = load(fullfile(root,'problem_5','data','linear_models.mat'), 'Ad','Bd','Ed','Cz');

%% Continuous-time observer canonical realizations (Problem 4.6)
% name, numerator gain, denominator, input (1 = F1, 2 = F2), output
blk = {'G11', Kid(1,1), [Tid(1) 1],                  1, 1;
       'G12', Kid(1,2), conv([Tid(1) 1],[Tid(3) 1]), 2, 1;
       'G21', Kid(2,1), conv([Tid(2) 1],[Tid(4) 1]), 1, 2;
       'G22', Kid(2,2), [Tid(2) 1],                  2, 2};
dgain = [0, 1/(1-gamma2), 1/(1-gamma1), 0];            % F3 -> G12 states, F4 -> G21 states
dcol  = [0, 1, 2, 0];

nb = size(blk,1);
nu = 2;  ny = 2;  nd = 2;
nxi = zeros(1,nb);
for i = 1:nb, nxi(i) = length(blk{i,3}) - 1; end
nx = sum(nxi);

Ac = zeros(nx);  Bc = zeros(nx,nu);  Cc = zeros(ny,nx);  Ec = zeros(nx,nd);
idx = [0 cumsum(nxi)];
for i = 1:nb
    den = blk{i,3}/blk{i,3}(1);                        % monic denominator
    b0  = blk{i,2}/blk{i,3}(1);                        % numerator (constant)
    n   = nxi(i);
    a   = den(2:end);
    Ai  = [-a(:), [eye(n-1); zeros(1,n-1)]];           % observer canonical form
    Bi  = [zeros(n-1,1); b0];                          % (num has no s-terms)
    Ci  = [1 zeros(1,n-1)];
    r   = idx(i)+1:idx(i+1);
    Ac(r,r)          = Ai;
    Bc(r,blk{i,4})   = Bi;
    Cc(blk{i,5},r)   = Ci;
    if dcol(i) > 0
        Ec(r,dcol(i)) = dgain(i)*Bi;
    end
end
C = Cc;

%% Discretization (ZOH) and process noise
Phi = expm([Ac Bc Ec; zeros(nu+nd, nx+nu+nd)]*Ts);
Ad  = Phi(1:nx, 1:nx);
Bd  = Phi(1:nx, nx+1:nx+nu);
Ed  = Phi(1:nx, nx+nu+1:end);

Gc        = Ec*diag(sigma);                            % SDE diffusion (6 x 2)
[~,Qsde]  = DiscretizeSDE(Ac, Gc, Ts);
Qpc       = Ed*Rdd*Ed';

% Checks: steady-state gains against Problem 4 and the Problem 5 model
Ku   = C*((eye(nx)-Ad)\Bd);
Kd   = C*((eye(nx)-Ad)\Ed);
Kd5  = P5.Cz*((eye(4)-P5.Ad)\P5.Ed);
fprintf('Gain check u -> y: max|C(I-Ad)^{-1}Bd - Kid| = %.1e\n', max(abs(Ku(:)-Kid(:))));
fprintf('Disturbance gain d -> y [cm/(cm^3/s)]:\n');
fprintf('   Problem 4 model: [%.4f %.4f; %.4f %.4f]\n', Kd');
fprintf('   Problem 5 model: [%.4f %.4f; %.4f %.4f]\n', Kd5');
Ob = zeros(ny*nx,nx);  Mk = C;
for k = 1:nx
    Ob((k-1)*ny+1:k*ny,:) = Mk;  Mk = Mk*Ad;
end
fprintf('Observability rank of (Ad,C): %d of %d (T1 and T2 appear in two blocks each)\n', ...
        rank(Ob), nx);

%% Disturbance models
mdl(1).name = 'Piecewise constant d_k';
mdl(1).tag  = 'piecewise';
mdl(1).G    = Ed;
mdl(1).Q    = Rdd;

mdl(2).name = 'SDE (Brownian motion)';
mdl(2).tag  = 'sde';
mdl(2).G    = eye(nx);
mdl(2).Q    = Qsde;

%% Time grid and inputs
tfinal = 40*60;                    % [s]
t      = 0:Ts:tfinal;
tmin   = t/60;
N      = length(t);

du = zeros(nu,N);                  % F1 = F2 = us = 300 cm^3/s

%% Initial condition (CHOICE)
% True plant at the steady state with both levels 2 cm above hs (reached by
% an earlier constant input offset); the filters start at the operating point.
du0   = Kid\[2; 2];
x0    = (eye(nx)-Ad)\(Bd*du0);
xhat0 = zeros(nx,1);

tconv = 5*60;
iT    = t <  tconv;
iS    = t >= tconv;

%% Simulation of both filters for both models
res = struct([]);
for m = 1:numel(mdl)
    G  = mdl(m).G;
    Q  = mdl(m).Q;

    kf = KalmanFilterDesign(Ad,Bd,C,G,Q,Rvv);
    P0 = diag((1.5*x0).^2) + kf.P;     % prior: 1.5 x the initial offset per state

    % noise realizations (same seed and order as sim_6_2_kalman_linear)
    rng(seed);
    v = chol(Rvv,'lower')*randn(ny,N);
    switch mdl(m).tag
        case 'piecewise'
            w = chol(Q,'lower')*randn(nd,N);
            D = ds + w;                                  % F3, F4 [cm^3/s]
        case 'sde'
            nsub  = 100;
            dt    = Ts/nsub;
            PhiF  = expm([Ac Bc; zeros(nu,nx+nu)]*dt);
            Af    = PhiF(1:nx,1:nx);
            Bf    = PhiF(1:nx,nx+1:end);
            dW    = sqrt(dt)*randn(2,nsub,N);
            tfine = (0:nsub*N)*dt;
            D     = sigma.*[zeros(2,1) cumsum(reshape(dW,2,[]),2)];   % sigma*omega [cm^3]
    end

    X    = zeros(nx,N);  Y   = zeros(ny,N);
    XfD  = zeros(nx,N);  XfS = zeros(nx,N);
    sdD  = zeros(ny,N);                      % sqrt(diag(C P_{k|k} C'))
    nisD = zeros(1,N);   nisS = zeros(1,N);
    Kerr = zeros(1,N);

    x   = x0;
    xpD = xhat0;  PpD = P0;
    xpS = xhat0;

    for k = 1:N
        X(:,k) = x;
        Y(:,k) = C*x + v(:,k);

        Kk      = (PpD*C')/(C*PpD*C' + Rvv);
        Kerr(k) = norm(Kk - kf.Kfx)/norm(kf.Kfx);

        [XfD(:,k),PfD,xpD,PpD,eD,ReD] = KalmanFilterDynamic(kf,xpD,PpD,Y(:,k),du(:,k));
        [XfS(:,k),xpS,eS]             = KalmanFilterStatic(kf,xpS,Y(:,k),du(:,k));

        sdD(:,k) = sqrt(diag(C*PfD*C'));
        nisD(k)  = eD'*(ReD\eD);
        nisS(k)  = eS'*(kf.Re\eS);

        % true (linear) plant
        switch mdl(m).tag
            case 'piecewise'
                x = Ad*x + Bd*du(:,k) + Ed*w(:,k);
            case 'sde'
                for j = 1:nsub
                    x = Af*x + Bf*du(:,k) + Gc*dW(:,j,k);
                end
        end
    end

    res(m).kf   = kf;
    res(m).H    = C*X;                       % true level deviation  [cm]
    res(m).Ym   = Y;                         % measured h1, h2       [cm]
    res(m).HfD  = C*XfD;                     % dynamic KF estimate   [cm]
    res(m).HfS  = C*XfS;                     % static  KF estimate   [cm]
    res(m).sdD  = sdD;
    res(m).sdS  = sqrt(diag(C*kf.Pf*C'));
    res(m).nisD = nisD;
    res(m).nisS = nisS;
    res(m).Kerr = Kerr;
    res(m).D    = D;
    if strcmp(mdl(m).tag,'sde')
        res(m).tD = tfine;
    end
end

%% Evaluation
rmse = @(E,idx) sqrt(mean(E(:,idx).^2,2));
for m = 1:numel(mdl)
    r  = res(m);
    eD = r.HfD - r.H;
    eS = r.HfS - r.H;
    eY = r.Ym  - r.H;

    fprintf('\n=== Problem 4 model, %s ===\n', mdl(m).name);
    fprintf('Stationary output gain C*Kfx (y_{k|k} = y_{k|k-1} + C Kfx e_k):\n');
    disp(C*r.kf.Kfx)
    fprintf('%-26s %8s %8s\n','RMSE [cm]','h1','h2');
    fprintf('%-26s %8.3f %8.3f\n','measurement, t >= 5 min', rmse(eY,iS));
    fprintf('%-26s %8.3f %8.3f\n','static KF,  t <  5 min',  rmse(eS,iT));
    fprintf('%-26s %8.3f %8.3f\n','dynamic KF, t <  5 min',  rmse(eD,iT));
    fprintf('%-26s %8.3f %8.3f\n','static KF,  t >= 5 min',  rmse(eS,iS));
    fprintf('%-26s %8.3f %8.3f\n','dynamic KF, t >= 5 min',  rmse(eD,iS));
    fprintf('%-26s %8.3f %8.3f\n','theory sqrt(diag(CPfC''))', r.sdS);
    fprintf('Mean NIS for t >= 5 min (expected %d): static %.2f, dynamic %.2f\n', ...
            ny, mean(r.nisS(iS)), mean(r.nisD(iS)));
    kc = find(r.Kerr < 0.01, 1);
    if ~isempty(kc)
        fprintf('Dynamic gain within 1%% of the stationary gain after %.1f min\n', tmin(kc));
    end
end

%% Plots
plotdir = fullfile(here,'plots');
if ~exist(plotdir,'dir'), mkdir(plotdir), end

cS = [0 0.447 0.741];
cD = [0.850 0.325 0.098];
cY = [0.6 0.6 0.6];
cF = [0.494 0.184 0.556];
hname = {'h_1','h_2'};
dname = {'F_3','F_4'};
figs  = {};
names = {};

% Figure 1/2: measured levels (row 1) and unknown disturbances (row 2)
for m = 1:numel(mdl)
    r = res(m);
    fig = figure('Name',['6.2 Problem 4 model, ' mdl(m).name],'Position',[100 100 1000 600]);
    for i = 1:ny
        subplot(2,2,i)
        plot(tmin, hs(i)+r.Ym(i,:),  '.',  'Color', cY, 'MarkerSize', 5); hold on
        plot(tmin, hs(i)+r.H(i,:),   'k-', 'LineWidth', 1.2);
        plot(tmin, hs(i)+r.HfS(i,:), '-',  'Color', cS, 'LineWidth', 1.0);
        plot(tmin, hs(i)+r.HfD(i,:), '--', 'Color', cD, 'LineWidth', 1.0);
        hold off; grid on; box on; xlim([0 tfinal/60])
        ylabel([hname{i} ' [cm]'])
        title([hname{i} ' - ' mdl(m).name ' (Problem 4 model)'])
        if i == 1
            legend('measurement','true','static KF','dynamic KF','Location','best')
        end

        subplot(2,2,2+i)
        if strcmp(mdl(m).tag,'sde')
            plot(r.tD/60, r.D(i,:), 'Color', cF, 'LineWidth', 0.8); hold on
            plot([0 tfinal/60], [0 0], 'k--', 'LineWidth', 1.0);
            hold off; grid on; box on; xlim([0 tfinal/60])
            ylabel(['\int(' dname{i} ' - d_{s,' num2str(i) '}) dt [cm^3]']); xlabel('time [min]')
            title([dname{i} ': integrated disturbance \sigma_{' dname{i} '}\omega_' num2str(i) '(t)'])
        else
            stairs(tmin, r.D(i,:), 'Color', cF, 'LineWidth', 0.8); hold on
            plot([0 tfinal/60], ds(i)*[1 1], 'k--', 'LineWidth', 1.0);
            hold off; grid on; box on; xlim([0 tfinal/60])
            ylabel([dname{i} ' [cm^3/s]']); xlabel('time [min]')
            title([dname{i} ' (held constant over each sample)'])
            legend('realization d_k','mean d_s','Location','best')
        end
    end
    figs{end+1}  = fig;
    names{end+1} = ['sim_6_2_p4_levels_' mdl(m).tag];
end

% Figure 3: estimation errors of h1, h2 with +-2 std bands
fig = figure('Name','6.2 Problem 4 model, estimation errors','Position',[100 100 1000 500]);
for m = 1:numel(mdl)
    r = res(m);
    for i = 1:ny
        subplot(ny,2,(i-1)*2+m)
        fill([tmin fliplr(tmin)], [2*r.sdD(i,:) -2*fliplr(r.sdD(i,:))], cD, ...
             'FaceAlpha', 0.15, 'EdgeColor', 'none'); hold on
        plot(tmin, r.HfS(i,:)-r.H(i,:), '-',  'Color', cS, 'LineWidth', 0.9);
        plot(tmin, r.HfD(i,:)-r.H(i,:), '--', 'Color', cD, 'LineWidth', 0.9);
        plot([0 tfinal/60], 2*r.sdS(i)*[1 1], 'k:', [0 tfinal/60], -2*r.sdS(i)*[1 1], 'k:');
        hold off; grid on; box on; xlim([0 tfinal/60])
        ylabel(['e_{' hname{i}(3) '} [cm]'])
        if i == 1
            title([mdl(m).name ' (Problem 4 model)'])
        end
        if i == 1 && m == 1
            legend('\pm2\sigma dynamic','static KF','dynamic KF','\pm2\sigma stationary','Location','best')
        end
        if i == ny, xlabel('time [min]'), end
    end
end
figs{end+1}  = fig;
names{end+1} = 'sim_6_2_p4_errors';

for f = 1:numel(figs)
    try, theme(figs{f},'light'), catch, end
    exportgraphics(figs{f}, fullfile(plotdir,[names{f} '.png']), 'Resolution', 150)
    exportgraphics(figs{f}, fullfile(plotdir,[names{f} '.pdf']), 'ContentType', 'vector')
end

%% Save the model and filter designs for later problems
datadir = fullfile(here,'data');
if ~exist(datadir,'dir'), mkdir(datadir), end
kfPiecewise = res(1).kf;
kfSDE       = res(2).kf;
save(fullfile(datadir,'kalman_6_2_problem4.mat'), ...
     'Ac','Bc','Cc','Ec','Gc','Ad','Bd','Ed','C','Qpc','Qsde', ...
     'kfPiecewise','kfSDE','Ts');
