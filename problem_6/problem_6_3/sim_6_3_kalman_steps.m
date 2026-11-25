%% Problem 6.3 - Kalman filters with step changes in the unknown disturbances
% The unknown disturbances F3, F4 are stochastic AND contain step changes:
%
%   d_k = ds + dstep_k + w_k,   w_k ~ N(0,Rdd)  (piecewise constant, default)
%
% with ds = [100; 80] cm^3/s and dstep a +20 cm^3/s step in F3 at 10 min and
% in F4 at 25 min. The pumps are held at us = [300; 300] cm^3/s (du = 0).
%
% Filters (static and dynamic versions of each):
%   plain      : the 6.2 design. Assumes zero-mean disturbances, so a step
%                is an unmodelled constant bias -> steady-state offset.
%   augmented  : integrating (random walk) disturbance model on F3, F4
%                  dbar_{k+1} = dbar_k + xi_k,  xi_k ~ N(0,Qxi)
%                  [x;dbar]_{k+1} = [A E; 0 I][x;dbar]_k + [B;0]u_k + [G 0;0 I][w;xi]_k
%                  y_k = [C 0][x;dbar]_k + v_k
%                The filter estimates the shifted mean of F3, F4 -> no offset.
%   input dist.: the alternative of Azam (2017), Sec. 3.5.1: random walk on
%                the inputs, [A B; 0 I] (static filter, offsets only in table).
%
% Both linear models are used: the linearized model from Problem 5
% (4 states) and the identified transfer function model from Problem 4
% (6-state realization, Problem4StateSpace). The plant is the same linear
% model as the filter (the nonlinear plant is Problem 6.4).
%
% Uses KalmanFilterDesign, KalmanFilterRun, Problem4StateSpace (../),
% DiscretizeSDE (problem_5).

clc
clear
close all

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..');
addpath(fullfile(here,'..'));                          % Kalman filter functions
addpath(fullfile(root,'problem_2','problem_2_1'));     % parameters
addpath(fullfile(root,'problem_5'));                   % DiscretizeSDE

QuadrupleTankParameters                                % rho, At, gamma1, gamma2, us, ds, Rdd, Rvv, sigma, seed
P5 = load(fullfile(root,'problem_5','data','linear_models.mat'));
P4 = load(fullfile(root,'problem_4','data','tf_identification.mat'), 'Kid', 'Tid');
Ts = P5.Ts;
hs = P5.hs;

%% Settings (CHOICE)
noiseModel = 'piecewise';          % 'piecewise' (Problem 2.2) or 'sde' (Problem 2.3)

tfinal = 40*60;                    % [s]
t      = 0:Ts:tfinal;
tmin   = t/60;
N      = length(t);
nu = 2;  ny = 2;  nd = 2;

du = zeros(nu,N);                  % F1 = F2 = 300 cm^3/s

tstep = [10 25]*60;                % [s] step times in F3 and F4
dstep = zeros(nd,N);               % step part of the disturbance deviation
dstep(1, t >= tstep(1)) = 20;      % F3: 100 -> 120 cm^3/s
dstep(2, t >= tstep(2)) = 20;      % F4:  80 -> 100 cm^3/s

sigXi     = 2;                     % [cm^3/s] std of the random walk per sample
Qxi       = sigXi^2*eye(nd);
sigXiTune = [0.5 2 8];             % tuning study
P0d       = 20^2*eye(nd);          % prior covariance of the disturbance states

% Windows for the steady-state offset (last 5 min before the next event)
win = {[5 10], [20 25], [35 40]};
winName = {'before steps (5-10 min)','after F3 step (20-25 min)','after both steps (35-40 min)'};

%% Linear models
m2h = 1./(rho*At);

mdl(1).name  = 'Problem 5 model (linearized)';
mdl(1).tag   = 'p5';
mdl(1).Ad = P5.Ad;  mdl(1).Bd = P5.Bd;  mdl(1).Ed = P5.Ed;  mdl(1).C = P5.C;
mdl(1).Ac = P5.A;   mdl(1).Bc = P5.B;   mdl(1).Ec = P5.E;
mdl(1).Gc = P5.G;   mdl(1).Qsde = P5.Qd;
mdl(1).H  = diag(m2h);                                % states -> levels h1..h4 [cm]
mdl(1).x0 = 2./m2h;                                   % 2 cm above hs in all tanks
mdl(1).P0 = diag((3./m2h).^2);

s4 = Problem4StateSpace(P4.Kid, P4.Tid, gamma1, gamma2, Ts);
mdl(2).name  = 'Problem 4 model (identified)';
mdl(2).tag   = 'p4';
mdl(2).Ad = s4.Ad;  mdl(2).Bd = s4.Bd;  mdl(2).Ed = s4.Ed;  mdl(2).C = s4.C;
mdl(2).Ac = s4.Ac;  mdl(2).Bc = s4.Bc;  mdl(2).Ec = s4.Ec;
mdl(2).Gc = s4.Ec*diag(sigma);
[~,mdl(2).Qsde] = DiscretizeSDE(s4.Ac, mdl(2).Gc, Ts);
mdl(2).H  = s4.C;                                     % states -> levels h1, h2 [cm]
mdl(2).x0 = (eye(size(s4.Ad))-s4.Ad)\(s4.Bd*(P4.Kid\[2; 2]));   % levels 2 cm above hs
mdl(2).P0 = [];                                       % set below from the stationary P

%% Simulation and filtering
res = struct([]);
for m = 1:numel(mdl)
    M  = mdl(m);
    nx = size(M.Ad,1);

    % ---------- noise realizations (same seed and order as Problem 6.2)
    rng(seed);
    v = chol(Rvv,'lower')*randn(ny,N);
    switch noiseModel
        case 'piecewise'
            G = M.Ed;  Q = Rdd;
            w = chol(Rdd,'lower')*randn(nd,N);
            D = ds + dstep + w;                       % realized F3, F4 [cm^3/s]
        case 'sde'
            G = eye(nx);  Q = M.Qsde;
            nsub = 100;  dt = Ts/nsub;
            Phi  = expm([M.Ac M.Bc M.Ec; zeros(nu+nd, nx+nu+nd)]*dt);
            Af   = Phi(1:nx,1:nx);
            Bf   = Phi(1:nx,nx+1:nx+nu);
            Ef   = Phi(1:nx,nx+nu+1:end);
            dW   = sqrt(dt)*randn(nd,nsub,N);
            D    = ds + dstep;                        % mean of F3, F4 (white noise on top)
    end

    % ---------- true (linear) plant
    X = zeros(nx,N);  Y = zeros(ny,N);
    x = M.x0;
    for k = 1:N
        X(:,k) = x;
        Y(:,k) = M.C*x + v(:,k);
        switch noiseModel
            case 'piecewise'
                x = M.Ad*x + M.Bd*du(:,k) + M.Ed*(dstep(:,k) + w(:,k));
            case 'sde'
                for j = 1:nsub
                    x = Af*x + Bf*du(:,k) + Ef*dstep(:,k) + M.Gc*dW(:,j,k);
                end
        end
    end

    % ---------- filter designs
    kfP = KalmanFilterDesign(M.Ad, M.Bd, M.C, G, Q, Rvv);         % plain (6.2)
    if isempty(M.P0)
        M.P0 = diag((1.5*M.x0).^2) + kfP.P;
    end

    aug = @(Adist,Qx) KalmanFilterDesign( ...
            [M.Ad Adist; zeros(nd,nx) eye(nd)], [M.Bd; zeros(nd,nu)], [M.C zeros(ny,nd)], ...
            blkdiag(G, eye(nd)), blkdiag(Q, Qx), Rvv);
    kfA = aug(M.Ed, Qxi);                                         % disturbance on F3, F4
    kfB = aug(M.Bd, Qxi);                                         % disturbance on inputs

    rP = KalmanFilterRun(kfP, Y, du, zeros(nx,1),    M.P0);
    rA = KalmanFilterRun(kfA, Y, du, zeros(nx+nd,1), blkdiag(M.P0, P0d));
    rB = KalmanFilterRun(kfB, Y, du, zeros(nx+nd,1), blkdiag(M.P0, P0d));

    % tuning study: augmented static filter for several Qxi
    tune = struct([]);
    for q = 1:numel(sigXiTune)
        kfT = aug(M.Ed, sigXiTune(q)^2*eye(nd));
        rT  = KalmanFilterRun(kfT, Y, du, zeros(nx+nd,1), blkdiag(M.P0, P0d));
        tune(q).H  = M.H*rT.XfS(1:nx,:);
        tune(q).dh = rT.XfS(nx+1:end,:);
    end

    % ---------- store (levels in cm, disturbance estimates in cm^3/s)
    nh = size(M.H,1);
    Haug = [M.H zeros(nh,nd)];
    r.H     = M.H*X;
    r.Ym    = Y;
    r.HPs   = M.H*rP.XfS;           r.HPd = M.H*rP.XfD;
    r.HAs   = Haug*rA.XfS;          r.HAd = Haug*rA.XfD;
    r.HBs   = Haug*rB.XfS;
    r.dAs   = rA.XfS(nx+1:end,:);   r.dAd = rA.XfD(nx+1:end,:);
    r.sdA   = zeros(nh,N);          r.sdP = zeros(nh,N);
    for k = 1:N
        r.sdA(:,k) = sqrt(diag(Haug*rA.PfD(:,:,k)*Haug'));
        r.sdP(:,k) = sqrt(diag(M.H*rP.PfD(:,:,k)*M.H'));
    end
    r.D     = D;
    r.tune  = tune;
    r.kfP   = kfP;  r.kfA = kfA;  r.kfB = kfB;
    r.nisP  = rP.nisS;  r.nisA = rA.nisS;
    res = [res r]; %#ok<AGROW>
    clear r
end

%% Steady-state offsets: mean estimation error in each window
for m = 1:numel(mdl)
    r  = res(m);
    nh = size(r.H,1);
    fprintf('\n=================== %s, %s disturbances ===================\n', mdl(m).name, noiseModel);
    filt = {'plain static',     r.HPs;
            'plain dynamic',    r.HPd;
            'augmented static', r.HAs;
            'augmented dynamic',r.HAd;
            'input-dist. static',r.HBs};
    for w = 1:numel(win)
        iw = tmin >= win{w}(1) & tmin <= win{w}(2);
        fprintf('\nMean estimation error [cm], %s\n', winName{w});
        fprintf('%-20s', 'filter');
        for i = 1:nh, fprintf('%9s', sprintf('h%d',i)); end
        fprintf('\n');
        for f = 1:size(filt,1)
            fprintf('%-20s', filt{f,1});
            fprintf('%9.3f', mean(filt{f,2}(:,iw) - r.H(:,iw), 2));
            fprintf('\n');
        end
        fprintf('%-20s', 'true F3, F4 mean'); fprintf('%9.1f', ds + mean(dstep(:,iw),2)); fprintf('\n');
        fprintf('%-20s', 'estimated (aug.)');  fprintf('%9.1f', ds + mean(r.dAs(:,iw),2));  fprintf('\n');
    end
    iS = tmin >= 5;
    fprintf('\nMean NIS for t >= 5 min (expected %d): plain %.2f, augmented %.2f\n', ...
            ny, mean(r.nisP(iS)), mean(r.nisA(iS)));
end

%% Plots
plotdir = fullfile(here,'plots');
if ~exist(plotdir,'dir'), mkdir(plotdir), end

cY  = [0.6 0.6 0.6];
cPs = [0 0.447 0.741];       cPd = [0.30 0.745 0.933];      % plain static / dynamic
cAs = [0.850 0.325 0.098];   cAd = [0.929 0.694 0.125];     % augmented static / dynamic
cF  = [0.75 0.65 0.85];
cT  = [0 0.6 0.3; 0.850 0.325 0.098; 0.494 0.184 0.556];    % tuning
dname = {'F_3','F_4'};
figs  = {};
names = {};

for m = 1:numel(mdl)
    r  = res(m);
    nh = size(r.H,1);

    % ---------- Figure 1: levels (row 1) and disturbances (row 2)
    fig = figure('Name',['6.3 ' mdl(m).name],'Position',[100 100 1100 650]);
    for i = 1:ny
        subplot(2,2,i)
        p0 = plot(tmin, hs(i)+r.Ym(i,:),  '.',  'Color', cY,  'MarkerSize', 5); hold on
        p1 = plot(tmin, hs(i)+r.H(i,:),   'k-', 'LineWidth', 1.3);
        p2 = plot(tmin, hs(i)+r.HPs(i,:), '-',  'Color', cPs, 'LineWidth', 1.0);
        p3 = plot(tmin, hs(i)+r.HPd(i,:), ':',  'Color', cPd, 'LineWidth', 1.2);
        p4 = plot(tmin, hs(i)+r.HAs(i,:), '-',  'Color', cAs, 'LineWidth', 1.0);
        p5 = plot(tmin, hs(i)+r.HAd(i,:), '--', 'Color', cAd, 'LineWidth', 1.0);
        yl = ylim;
        for s = 1:numel(tstep)
            plot(tstep(s)/60*[1 1], yl, 'k:', 'LineWidth', 0.8);
        end
        ylim(yl); hold off; grid on; box on; xlim([0 tfinal/60])
        ylabel(sprintf('h_%d [cm]', i))
        title(sprintf('h_%d - %s', i, mdl(m).name))
        if i == 1
            legend([p0 p1 p2 p3 p4 p5], {'measurement','true','plain static','plain dynamic', ...
                   'augmented static','augmented dynamic'}, 'Location','best')
        end

        subplot(2,2,2+i)
        q = [];
        if strcmp(noiseModel,'piecewise')
            q(end+1) = stairs(tmin, r.D(i,:), 'Color', cF, 'LineWidth', 0.6); hold on
        end
        q(end+1) = stairs(tmin, ds(i)+dstep(i,:), 'k-', 'LineWidth', 1.3); hold on
        q(end+1) = plot(tmin, ds(i)+r.dAs(i,:), '-',  'Color', cAs, 'LineWidth', 1.1);
        q(end+1) = plot(tmin, ds(i)+r.dAd(i,:), '--', 'Color', cAd, 'LineWidth', 1.1);
        hold off; grid on; box on; xlim([0 tfinal/60])
        ylabel([dname{i} ' [cm^3/s]']); xlabel('time [min]')
        title([dname{i} ': true mean and augmented-filter estimate'])
        if strcmp(noiseModel,'piecewise')
            lg = {'realization','true mean d_s + step','estimate (static)','estimate (dynamic)'};
        else
            lg = {'true mean d_s + step','estimate (static)','estimate (dynamic)'};
        end
        legend(q, lg, 'Location','best')
    end
    figs{end+1}  = fig;
    names{end+1} = ['sim_6_3_' mdl(m).tag '_levels_' noiseModel];

    % ---------- Figure 2: estimation errors of all physical levels
    fig = figure('Name',['6.3 errors, ' mdl(m).name],'Position',[100 100 1000 180*nh+120]);
    for i = 1:nh
        subplot(nh,1,i)
        yl = [-1 1]*max(1.5, 1.1*max(abs(r.HPs(i,tmin>=5)-r.H(i,tmin>=5))));
        band = min(2*r.sdA(i,:), yl(2));                 % clip the initial (large) band
        pb = fill([tmin fliplr(tmin)], [band -fliplr(band)], cAs, ...
                  'FaceAlpha', 0.12, 'EdgeColor', 'none'); hold on
        e1 = plot(tmin, r.HPs(i,:)-r.H(i,:), '-',  'Color', cPs, 'LineWidth', 1.0);
        e2 = plot(tmin, r.HPd(i,:)-r.H(i,:), ':',  'Color', cPd, 'LineWidth', 1.2);
        e3 = plot(tmin, r.HAs(i,:)-r.H(i,:), '-',  'Color', cAs, 'LineWidth', 1.0);
        e4 = plot(tmin, r.HAd(i,:)-r.H(i,:), '--', 'Color', cAd, 'LineWidth', 1.0);
        plot([0 tfinal/60], [0 0], 'k-', 'LineWidth', 0.5);
        for s = 1:numel(tstep)
            plot(tstep(s)/60*[1 1], yl, 'k:', 'LineWidth', 0.8);
        end
        ylim(yl); hold off; grid on; box on; xlim([0 tfinal/60])
        ylabel(sprintf('e_{h%d} [cm]', i))
        if i == 1
            title(['Estimation error, ' mdl(m).name])
            legend([pb e1 e2 e3 e4], {'\pm2\sigma augmented (dynamic)','plain static','plain dynamic', ...
                   'augmented static','augmented dynamic'}, 'Location','best')
        end
        if i > ny
            title(sprintf('h_%d (not measured)', i))
        end
        if i == nh, xlabel('time [min]'), end
    end
    figs{end+1}  = fig;
    names{end+1} = ['sim_6_3_' mdl(m).tag '_errors_' noiseModel];

    % ---------- Figure 3: tuning of the random-walk variance Qxi
    fig = figure('Name',['6.3 Qxi tuning, ' mdl(m).name],'Position',[100 100 1000 600]);
    for i = 1:nd
        subplot(2,2,i)
        q = stairs(tmin, ds(i)+dstep(i,:), 'k-', 'LineWidth', 1.3); hold on
        lg = {'true mean'};
        for qq = 1:numel(sigXiTune)
            q(end+1) = plot(tmin, ds(i)+r.tune(qq).dh(i,:), 'Color', cT(qq,:), 'LineWidth', 1.0); %#ok<SAGROW>
            lg{end+1} = sprintf('\\sigma_\\xi = %g', sigXiTune(qq)); %#ok<SAGROW>
        end
        hold off; grid on; box on; xlim([0 tfinal/60])
        ylabel([dname{i} ' [cm^3/s]'])
        title([dname{i} ' estimate (augmented static KF)'])
        if i == 1, legend(q, lg, 'Location','best'), end

        subplot(2,2,2+i)
        for qq = 1:numel(sigXiTune)
            plot(tmin, r.tune(qq).H(i,:)-r.H(i,:), 'Color', cT(qq,:), 'LineWidth', 0.9); hold on
        end
        plot([0 tfinal/60], [0 0], 'k-', 'LineWidth', 0.5);
        hold off; grid on; box on; xlim([0 tfinal/60]); ylim([-1.5 1.5])
        ylabel(sprintf('e_{h%d} [cm]', i)); xlabel('time [min]')
        title(sprintf('h_%d estimation error', i))
    end
    figs{end+1}  = fig;
    names{end+1} = ['sim_6_3_' mdl(m).tag '_tuning_' noiseModel];
end

for f = 1:numel(figs)
    try, theme(figs{f},'light'), catch, end
    exportgraphics(figs{f}, fullfile(plotdir,[names{f} '.png']), 'Resolution', 150)
    exportgraphics(figs{f}, fullfile(plotdir,[names{f} '.pdf']), 'ContentType', 'vector')
end

%% Save the offset-free filter designs (augmented, disturbance on F3/F4) for MPC
datadir = fullfile(here,'data');
if ~exist(datadir,'dir'), mkdir(datadir), end
kfAugP5 = res(1).kfA;
kfAugP4 = res(2).kfA;
save(fullfile(datadir,['kalman_6_3_' noiseModel '.mat']), 'kfAugP5', 'kfAugP4', 'Qxi', 'Ts', 'noiseModel');
