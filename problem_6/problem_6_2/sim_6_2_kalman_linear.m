%% Problem 6.2 - Static and dynamic Kalman filters, Problem 5 models, no disturbance steps
% Discrete-time linear models from Problem 5.4 (data/linear_models.mat),
% deviation variables around the operating point (xs, us, ds):
%
%   x_{k+1} = Ad x_k + Bd u_k + G w_k,   w_k ~ N(0,Q)
%   y_k     = C  x_k + v_k,              v_k ~ N(0,Rvv)
%
% Two disturbance models are compared (the "unknown disturbances" F3, F4):
%   1) piecewise constant d_k (Problem 2.2):  G = Ed,   Q = Rdd
%      -> F3, F4 are white noise around ds, held constant over each sample
%   2) SDE / Brownian motion  (Problem 2.3):  G = I,    Q = Qd (Van Loan)
%      -> the integrated effect of dw over one sample
% The disturbances are zero-mean stochastic variables WITHOUT step changes,
% so the model the filters are designed on is the true data-generating model.
%
% The plant is the same linear model (the nonlinear plant is Problem 6.4).
% F1 = F2 = 300 cm^3/s (operating point); F3, F4 fluctuate around
% ds = [100; 80] cm^3/s with no change in their mean.
% The filter is started from a wrong initial estimate so that the transient
% of the dynamic (time-varying) filter can be compared with the static
% (stationary) filter, cf. Azam (2017), Sec. 3.2-3.4.
%
% Uses KalmanFilterDesign, KalmanFilterDynamic, KalmanFilterStatic (../).

clc
clear
close all

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..');
addpath(fullfile(here,'..'));                          % Kalman filter functions
addpath(fullfile(root,'problem_2','problem_2_1'));     % parameters

QuadrupleTankParameters                                % rho, At, seed
load(fullfile(root,'problem_5','data','linear_models.mat'), ...
     'Ts','Ad','Bd','Ed','C','Qd','Rdd','Rvv','us','ds','hs');
Sc = load(fullfile(root,'problem_5','data','linear_models.mat'),'A','B','G');  % continuous-time

nx = size(Ad,1);
nu = size(Bd,2);
ny = size(C,1);
m2h = 1./(rho*At);                 % mass [g] -> level [cm]

%% Disturbance models (CHOICE: the two stochastic models of Problem 5)
mdl(1).name = 'Piecewise constant d_k';
mdl(1).tag  = 'piecewise';
mdl(1).G    = Ed;
mdl(1).Q    = Rdd;

mdl(2).name = 'SDE (Brownian motion)';
mdl(2).tag  = 'sde';
mdl(2).G    = eye(nx);
mdl(2).Q    = Qd;

%% Time grid and inputs (deviation variables)
tfinal = 40*60;                    % [s]
t      = 0:Ts:tfinal;
tmin   = t/60;
N      = length(t);

% Inputs held at the operating point, F1 = F2 = us = 300 cm^3/s
% (the system is excited only by the stochastic disturbances F3, F4)
du = zeros(nu,N);

%% Initial condition (CHOICE)
% True state 2 cm above the operating point in all tanks; the filters start
% at the operating point with a prior std of 3 cm.
x0    = (2*ones(nx,1))./m2h;
xhat0 = zeros(nx,1);
P0    = diag((3./m2h).^2);

tconv = 5*60;                      % [s] end of the initial transient
iT    = t <  tconv;
iS    = t >= tconv;

%% Simulation of both filters for both models
res = struct([]);
for m = 1:numel(mdl)
    G = mdl(m).G;
    Q = mdl(m).Q;
    nw = size(Q,1);

    kf = KalmanFilterDesign(Ad,Bd,C,G,Q,Rvv);

    % noise realizations (same seed for both models)
    rng(seed);
    v = chol(Rvv,'lower')*randn(ny,N);
    switch mdl(m).tag
        case 'piecewise'
            % d_k = ds + w_k, w_k ~ N(0,Rdd), held constant over each sample
            w = chol(Q,'lower')*randn(nw,N);
            D = ds + w;                                  % F3, F4 [cm^3/s]
        case 'sde'
            % Brownian increments on a fine grid inside each sample; the plant
            % dX = (A X + B U) dt + G dw is integrated with exact drift and
            % Euler-Maruyama noise. Its sampled noise has covariance -> Qd.
            nsub = 100;
            dt   = Ts/nsub;
            Phi  = expm([Sc.A Sc.B; zeros(nu,nx+nu)]*dt);
            Af   = Phi(1:nx,1:nx);
            Bf   = Phi(1:nx,nx+1:end);
            dW   = sqrt(dt)*randn(2,nsub,N);
            % F3, F4 are white noise around ds; the Brownian motion is the
            % integrated disturbance  int (F - ds) dt = sigma*omega(t)  [cm^3]
            tfine = (0:nsub*N)*dt;
            D     = sigma.*[zeros(2,1) cumsum(reshape(dW,2,[]),2)];
    end

    X    = zeros(nx,N);  Y   = zeros(ny,N);
    XfD  = zeros(nx,N);  XfS = zeros(nx,N);
    sdD  = zeros(nx,N);                      % sqrt(diag(P_{k|k})), dynamic
    nisD = zeros(1,N);   nisS = zeros(1,N);  % normalized innovation squared
    Kerr = zeros(1,N);                       % ||K_k - K|| / ||K||

    x   = x0;
    xpD = xhat0;  PpD = P0;                  % dynamic:  x_{0|-1}, P_{0|-1}
    xpS = xhat0;                             % static:   x_{0|-1}

    for k = 1:N
        X(:,k) = x;
        Y(:,k) = C*x + v(:,k);

        Kk      = (PpD*C')/(C*PpD*C' + Rvv);
        Kerr(k) = norm(Kk - kf.Kfx)/norm(kf.Kfx);

        [XfD(:,k),PfD,xpD,PpD,eD,ReD] = KalmanFilterDynamic(kf,xpD,PpD,Y(:,k),du(:,k));
        [XfS(:,k),xpS,eS]             = KalmanFilterStatic(kf,xpS,Y(:,k),du(:,k));

        sdD(:,k) = sqrt(diag(PfD));
        nisD(k)  = eD'*(ReD\eD);
        nisS(k)  = eS'*(kf.Re\eS);

        % true (linear) plant
        switch mdl(m).tag
            case 'piecewise'
                x = Ad*x + Bd*du(:,k) + G*w(:,k);
            case 'sde'
                for j = 1:nsub
                    x = Af*x + Bf*du(:,k) + Sc.G*dW(:,j,k);
                end
        end
    end

    res(m).kf   = kf;
    res(m).H    = X  .*m2h;                  % true level deviation  [cm]
    res(m).Ym   = Y;                         % measured h1, h2       [cm]
    res(m).HfD  = XfD.*m2h;                  % dynamic KF estimate   [cm]
    res(m).HfS  = XfS.*m2h;                  % static  KF estimate   [cm]
    res(m).sdD  = sdD.*m2h;                  % dynamic KF std        [cm]
    res(m).sdS  = sqrt(diag(kf.Pf)).*m2h;    % stationary std        [cm]
    res(m).nisD = nisD;
    res(m).nisS = nisS;
    res(m).Kerr = Kerr;
    res(m).D    = D;                         % piecewise: F3, F4 [cm^3/s]; SDE: sigma*omega [cm^3]
    if strcmp(mdl(m).tag,'sde')
        res(m).tD = tfine;
    end
end

%% Evaluation: RMSE of the filtered levels
rmse = @(E,idx) sqrt(mean(E(:,idx).^2,2));
for m = 1:numel(mdl)
    r = res(m);
    eD = r.HfD - r.H;
    eS = r.HfS - r.H;
    eY = r.Ym  - r.H(1:ny,:);

    fprintf('\n=== %s ===\n', mdl(m).name);
    fprintf('Stationary gain Kfx (x_{k|k} = x_{k|k-1} + Kfx e_k), level units [cm/cm]:\n');
    disp(diag(m2h)*r.kf.Kfx)
    fprintf('%-26s %8s %8s %8s %8s\n','RMSE [cm]','h1','h2','h3','h4');
    fprintf('%-26s %8.3f %8.3f %8s %8s\n','measurement, t >= 5 min', rmse(eY,iS), '-', '-');
    fprintf('%-26s %8.3f %8.3f %8.3f %8.3f\n','static KF,  t <  5 min',  rmse(eS,iT));
    fprintf('%-26s %8.3f %8.3f %8.3f %8.3f\n','dynamic KF, t <  5 min',  rmse(eD,iT));
    fprintf('%-26s %8.3f %8.3f %8.3f %8.3f\n','static KF,  t >= 5 min',  rmse(eS,iS));
    fprintf('%-26s %8.3f %8.3f %8.3f %8.3f\n','dynamic KF, t >= 5 min',  rmse(eD,iS));
    fprintf('%-26s %8.3f %8.3f %8.3f %8.3f\n','theory sqrt(diag(Pf))',   r.sdS);
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

cS = [0 0.447 0.741];              % static
cD = [0.850 0.325 0.098];          % dynamic
cY = [0.6 0.6 0.6];                % measurement
hname = {'h_1','h_2','h_3','h_4'};
figs  = {};
names = {};

% Figure 1/2: measured levels (row 1) and unknown disturbances (row 2)
dname = {'F_3','F_4'};
for m = 1:numel(mdl)
    r = res(m);
    fig = figure('Name',['6.2 ' mdl(m).name],'Position',[100 100 1000 600]);
    for i = 1:ny
        subplot(2,2,i)
        plot(tmin, hs(i)+r.Ym(i,:),  '.',  'Color', cY, 'MarkerSize', 5); hold on
        plot(tmin, hs(i)+r.H(i,:),   'k-', 'LineWidth', 1.2);
        plot(tmin, hs(i)+r.HfS(i,:), '-',  'Color', cS, 'LineWidth', 1.0);
        plot(tmin, hs(i)+r.HfD(i,:), '--', 'Color', cD, 'LineWidth', 1.0);
        hold off; grid on; box on; xlim([0 tfinal/60])
        ylabel([hname{i} ' [cm]'])
        title(sprintf('%s - %s', hname{i}, mdl(m).name))
        if i == 1
            legend('measurement','true','static KF','dynamic KF','Location','best')
        end

        subplot(2,2,2+i)
        if strcmp(mdl(m).tag,'sde')
            % Brownian motion: integrated deviation of F3/F4 from the mean
            plot(r.tD/60, r.D(i,:), 'Color', [0.494 0.184 0.556], 'LineWidth', 0.8); hold on
            plot([0 tfinal/60], [0 0], 'k--', 'LineWidth', 1.0);
            hold off; grid on; box on; xlim([0 tfinal/60])
            ylabel(['\int(' dname{i} ' - d_{s,' num2str(i) '}) dt [cm^3]']); xlabel('time [min]')
            title([dname{i} ': integrated disturbance \sigma_{' dname{i} '}\omega_' num2str(i) '(t)'])
        else
            stairs(tmin, r.D(i,:), 'Color', [0.494 0.184 0.556], 'LineWidth', 0.8); hold on
            plot([0 tfinal/60], ds(i)*[1 1], 'k--', 'LineWidth', 1.0);
            hold off; grid on; box on; xlim([0 tfinal/60])
            ylabel([dname{i} ' [cm^3/s]']); xlabel('time [min]')
            title(sprintf('%s (held constant over each sample)', dname{i}))
            legend('realization d_k','mean d_s','Location','best')
        end
    end

    figs{end+1}  = fig;
    names{end+1} = ['sim_6_2_levels_' mdl(m).tag];
end

% Figure 3: estimation errors with +-2 std bands (dynamic) and stationary bands
fig = figure('Name','6.2 estimation errors','Position',[100 100 1000 800]);
for m = 1:numel(mdl)
    r = res(m);
    for i = 1:nx
        subplot(nx,2,(i-1)*2+m)
        fill([tmin fliplr(tmin)], [2*r.sdD(i,:) -2*fliplr(r.sdD(i,:))], cD, ...
             'FaceAlpha', 0.15, 'EdgeColor', 'none'); hold on
        plot(tmin, r.HfS(i,:)-r.H(i,:), '-',  'Color', cS, 'LineWidth', 0.9);
        plot(tmin, r.HfD(i,:)-r.H(i,:), '--', 'Color', cD, 'LineWidth', 0.9);
        plot([0 tfinal/60], 2*r.sdS(i)*[1 1], 'k:', [0 tfinal/60], -2*r.sdS(i)*[1 1], 'k:');
        hold off; grid on; box on; xlim([0 tfinal/60])
        ylabel(['e_{' hname{i}(3) '} [cm]'])
        if i == 1
            title(mdl(m).name)
        end
        if i == 1 && m == 1
            legend('\pm2\sigma dynamic','static KF','dynamic KF','\pm2\sigma stationary','Location','best')
        end
        if i == nx, xlabel('time [min]'), end
    end
end
figs{end+1}  = fig;
names{end+1} = 'sim_6_2_errors';

for f = 1:numel(figs)
    try, theme(figs{f},'light'), catch, end
    exportgraphics(figs{f}, fullfile(plotdir,[names{f} '.png']), 'Resolution', 150)
    exportgraphics(figs{f}, fullfile(plotdir,[names{f} '.pdf']), 'ContentType', 'vector')
end

%% Save the filter designs for later problems
datadir = fullfile(here,'data');
if ~exist(datadir,'dir'), mkdir(datadir), end
kfPiecewise = res(1).kf;
kfSDE       = res(2).kf;
save(fullfile(datadir,'kalman_6_2.mat'), 'kfPiecewise', 'kfSDE', 'mdl', 'Ts');
