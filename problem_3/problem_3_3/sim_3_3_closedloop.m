%% Problem 3.2 / 3.3 - decentralized P, PI and PID control of the four tank system
% Closed-loop simulation of three controllers on the three models of Problem 2.
%   rows:    P, PI, PID
%   columns: 2.1 deterministic, 2.2 piecewise constant stochastic, 2.3 SDE
% Pairing from 3.1: F1 <-> h1, F2 <-> h2.
% Tuning: SIMC rules (Skogestad, 2003) on the first-order direct channels.

clc
clear
close all

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','..','problem_2','problem_2_1'));
addpath(fullfile(here,'..','..','problem_2','problem_2_4'));
addpath(fullfile(here,'..','problem_3_2'));

QuadrupleTankParameters   % p, us, ds, xs, zs, hs, Rdd, Rvv, sigma, seed, ...

%% Time grid
Ts     = 10;        % [s]
tfinal = 40*60;     % [s]
t      = 0:Ts:tfinal;
tmin   = t/60;
N      = length(t);

nx = 4; nu = 2; ny = 2;
nsub = 10;          % Euler-Maruyama sub-steps per interval (SDE model)

%% Scenario  (CHOICE)
% setpoint steps: +5 cm in h1 at 2 min, +5 cm in h2 at 15 min
r = repmat(zs,1,N);
r(1, t >= 2*60)  = zs(1) + 5;
r(2, t >= 15*60) = zs(2) + 5;

% disturbances held at their nominal values
dnom = repmat(ds,1,N);

% pump limits
umin = [0; 0];
umax = [500; 500];

%% 3.2 SIMC tuning on the direct channels, K/(T s + 1)
T    = (At./ap).*sqrt(2*hs/g);               % tank time constants [s]
Kp   = [gamma1*T(1)/A1; gamma2*T(2)/A2];     % direct-channel gains [cm per cm^3/s]
Tp   = T(1:2);                               % process time constants of the paired loops
tauc = 0.5*Tp;                               % desired closed-loop time constants  (CHOICE)

Kc   = Tp./(Kp.*tauc);                       % SIMC proportional gain
tauI = min(Tp, 4*tauc);                      % SIMC integral time
tauD = tauI/4;                               % derivative time for the PID case   (CHOICE)
Nf   = [10; 10];                             % derivative filter ratio

ctrlNames = {'P','PI','PID'};
tunings   = {struct('Kc',Kc,'tauI',[Inf;Inf],'tauD',[0;0]), ...
             struct('Kc',Kc,'tauI',tauI,     'tauD',[0;0]), ...
             struct('Kc',Kc,'tauI',tauI,     'tauD',tauD)};

fprintf('SIMC tuning (loop 1: F1->h1, loop 2: F2->h2)\n')
fprintf('  K    = [%.4f %.4f] cm/(cm^3/s)\n', Kp)
fprintf('  T    = [%.1f %.1f] s\n', Tp)
fprintf('  tauc = [%.1f %.1f] s\n', tauc)
fprintf('  Kc   = [%.2f %.2f] (cm^3/s)/cm\n', Kc)
fprintf('  tauI = [%.1f %.1f] s\n', tauI)
fprintf('  tauD = [%.1f %.1f] s (PID only)\n\n', tauD)

%% Noise realizations, shared by all controllers so the columns are comparable
[W,Tw,~] = StdWienerProcess(tfinal,(N-1)*nsub,2,1,seed);   % 2.3 Wiener path (seeds rng)
w = chol(Rdd,'lower')*randn(2,N);                           % 2.2 disturbance noise
v = chol(Rvv,'lower')*randn(ny,N);                          % measurement noise (2.2, 2.3)
psde = [p; sigma];

%% Closed-loop simulations: 3 controllers x 3 models
modelNames = {'2.1 deterministic','2.2 piecewise stochastic','2.3 SDE'};
modelKeys  = {'deterministic','piecewise','sde'};
res = cell(3,3);

for ic = 1:3
    for im = 1:3
        ctrl      = tunings{ic};
        ctrl.N    = Nf;
        ctrl.ubar = us;
        ctrl.umin = umin;
        ctrl.umax = umax;

        res{ic,im} = closedloop(modelKeys{im},ctrl,r,dnom,w,v,W,Tw,nsub,t,Ts,xs,p,psde);
    end
end

%% Performance: integral absolute error of the true levels and total input movement
fprintf('%-6s %-26s %10s %10s %12s\n','ctrl','model','IAE h1','IAE h2','sum|du|')
for ic = 1:3
    for im = 1:3
        R  = res{ic,im};
        iae = Ts*sum(abs(r - R.h(1:2,:)),2)/60;             % [cm*min]
        mov = sum(abs(diff(R.u,1,2)),2);                      % [cm^3/s]
        fprintf('%-6s %-26s %10.1f %10.1f %12.0f\n', ...
            ctrlNames{ic},modelNames{im},iae(1),iae(2),sum(mov))
    end
end

%% Plot: outputs (rows = controllers, columns = models)
figure('Name','3.3 outputs','Position',[100 100 1400 900]);
for ic = 1:3
    for im = 1:3
        R = res{ic,im};
        subplot(3,3,(ic-1)*3+im)
        plot(tmin,R.h(1,:),'LineWidth',1.2)
        hold on
        plot(tmin,R.h(2,:),'LineWidth',1.2)
        stairs(tmin,r(1,:),'k--','LineWidth',0.8)
        stairs(tmin,r(2,:),'k--','LineWidth',0.8)
        ylabel('[cm]')
        if ic == 3, xlabel('time [min]'), end
        title(sprintf('%s on %s',ctrlNames{ic},modelNames{im}))
        if ic == 1 && im == 1
            legend('h_1','h_2','setpoints','Location','best')
        end
        grid on
    end
end

%% Plot: inputs
figure('Name','3.3 inputs','Position',[100 100 1400 900]);
for ic = 1:3
    for im = 1:3
        R = res{ic,im};
        subplot(3,3,(ic-1)*3+im)
        stairs(tmin,R.u(1,:),'LineWidth',1.2)
        hold on
        stairs(tmin,R.u(2,:),'LineWidth',1.2)
        ylabel('[cm^3/s]')
        if ic == 3, xlabel('time [min]'), end
        title(sprintf('%s on %s',ctrlNames{ic},modelNames{im}))
        if ic == 1 && im == 1
            legend('F_1','F_2','Location','best')
        end
        grid on
    end
end

%% Save figures to problem_3/plots
plotdir = fullfile(here,'..','plots');
if ~exist(plotdir,'dir')
    mkdir(plotdir)
end
names = {'sim_3_3_outputs','sim_3_3_inputs'};
for i = 1:2
    f = figure(i);
    try, theme(f,'light'), catch, end
    exportgraphics(f,fullfile(plotdir,[names{i} '.png']),'Resolution',150)
    try
        exportgraphics(f,fullfile(plotdir,[names{i} '.pdf']),'ContentType','vector')
    catch ME
        warning('PDF export of %s failed: %s',names{i},ME.message)
    end
end

%% Closed-loop simulation of one controller on one model
function out = closedloop(model,ctrl,r,dnom,w,v,W,Tw,nsub,t,Ts,xs,p,psde)

    N  = length(t);
    nx = length(xs);
    At  = p(5:8);
    rho = p(12);

    switch model
        case 'deterministic'
            d = dnom;       noisy = false;
        case 'piecewise'
            d = dnom + w;   noisy = true;    % d_k held constant over each interval
        case 'sde'
            d = dnom;       noisy = true;    % fluctuations enter through the diffusion
    end

    x = zeros(nx,N);
    y = zeros(2,N);
    u = zeros(2,N);
    x(:,1) = xs;

    for k = 1:N-1
        % measurement
        y(:,k) = QuadrupleTankSensor(x(:,k),p);
        if noisy
            y(:,k) = y(:,k) + v(:,k);
        end

        % controller
        [u(:,k),ctrl] = PIDController(ctrl,r(:,k),y(:,k),Ts);

        % plant over one sampling interval
        switch model
            case {'deterministic','piecewise'}
                [~,Xk] = ode15s(@QuadrupleTankProcess, ...
                                [t(k) t(k+1)], ...
                                x(:,k), ...
                                [], ...
                                u(:,k), d(:,k), p);
                x(:,k+1) = Xk(end,:)';
            case 'sde'
                idx = (k-1)*nsub + (1:nsub+1);
                Xk  = SDEsolverExplicitExplicit(@QuadrupleTankProcess, ...
                                                @QuadrupleTankDiffusion, ...
                                                Tw(idx), x(:,k), W(:,idx), ...
                                                u(:,k), d(:,k), psde);
                x(:,k+1) = Xk(:,end);
        end
    end

    y(:,N) = QuadrupleTankSensor(x(:,N),p);
    if noisy
        y(:,N) = y(:,N) + v(:,N);
    end
    u(:,N) = u(:,N-1);

    out.x = x;
    out.h = x./(rho*At);   % true levels
    out.y = y;
    out.u = u;
    out.d = d;
end
