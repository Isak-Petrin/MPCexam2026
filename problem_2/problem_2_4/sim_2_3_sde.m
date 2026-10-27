% Week 3 "Slide recreations.m": StdWienerProcess (03B slide 6),

clc
clear
close all

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','problem_2_1'));
addpath(fullfile(here,'..','problem_2_3'));

QuadrupleTankParameters   % p, us, ds, xs, Rvv, sigma, seed, ...

%% Time grid
Ts     = 10;        % [s] sampling time (measurements, input changes)
tfinal = 20*60;     % [s]
t      = 0:Ts:tfinal;
tmin   = t/60;
N      = length(t);

nsub = 10;          % Euler-Maruyama sub-steps per sampling interval

nx = 4;  % Number of states
ny = 2;  % Number of measurements
nz = 2;  % Number of controlled outputs

nw = 2;  % Two independent Wiener processes: B1 and B2

%% Piecewise constant input sequences (same as sim_2_1_deterministic)
F1 = [330 235 295 303 242 292 300 243 325 324 ...
      330 357 278 274 202 301 320 334 276 308];

F2 = [285 301 295 313 314 316 289 278 302 293 ...
      304 296 307 277 319 305 300 306 307 297];

U = [F1; F2];
u = [repelem(U,1,60/Ts) U(:,end)];

% No deterministic disturbance flows in the SDE drift.
% Stochastic forcing is introduced through Sigma*dW.
d = zeros(2,N);

% In case we want deterministic step incrase later:
% Mean disturbance with a step in F3 after 10 min
% d = repmat(ds,1,N);
% d(1, t >= 10*60) = ds(1) + 50;

%% Wiener process for the whole horizon (nw channels, one realization)
[W,Tw,~] = StdWienerProcess(tfinal,(N-1)*nsub,nw,1,seed);   % W: nw x ((N-1)*nsub+1)

%% Measurement noise
Lr = chol(Rvv,'lower');
v  = Lr*randn(ny,N);

%% Simulation: Euler-Maruyama over each sampling interval
psde = [p; sigma];   % psde(13:14) = sigma_F3, sigma_F4 for QuadrupleTankDiffusion

x = zeros(nx,N);
y = zeros(ny,N);
z = zeros(nz,N);
%  x(:,1) = xs; this used the initial state from models 1 and 2 - but xs is not the same
% for the SDE model as disturbances are removed from deterministic drift and enter through diffusion


%computing new steady state
% Use nominal pump flows and zero deterministic disturbances
dSDE = zeros(2,1);

u0 = u(:,1);


qsSDE = zeros(4,1);
qsSDE(3) = (1-gamma2)*u0(2);
qsSDE(4) = (1-gamma1)*u0(1);
qsSDE(1) = gamma1*u0(1) + qsSDE(3);
qsSDE(2) = gamma2*u0(2) + qsSDE(4);

hsSDE = (qsSDE./ap).^2/(2*g);
xsSDE = rho*At.*hsSDE;

x(:,1) = xsSDE;

for k = 1:N-1
    y(:,k) = QuadrupleTankSensor(x(:,k),p) + v(:,k);   % y = g(x,p) + v
    z(:,k) = QuadrupleTankOutput(x(:,k),p);

    idx = (k-1)*nsub + (1:nsub+1);                     % sub-grid of interval k

    Xk = SDEsolverExplicitExplicit(@QuadrupleTankProcess, ...
                                   @QuadrupleTankDiffusion, ...
                                   Tw(idx), ...
                                   x(:,k), ...
                                   W(:,idx), ...
                                   u(:,k), d(:,k), psde);

    x(:,k+1) = Xk(:,end);
end

y(:,N) = QuadrupleTankSensor(x(:,N),p) + v(:,N);
z(:,N) = QuadrupleTankOutput(x(:,N),p);

h = x./(rho*At);   % true levels [cm]


%% Plots: Brownian motions and measured vs. true tank levels
figure('Name','2.3 SDE: Brownian motions and tank levels');

% Top: Two independent Brownian motions
subplot(2,1,1)
plot(Tw/60, W(1,:), 'LineWidth',1.2)
hold on
plot(Tw/60, W(2,:), 'LineWidth',1.2)
xlabel('Time [min]')
ylabel('W(t)')
title('Standard Wiener Processes')
legend('B_1(t)','B_2(t)','Location','best')
grid on

% Bottom: Measured and true levels for tanks 1 and 2
subplot(2,1,2)
plot(tmin, y(1,:), '-', 'LineWidth',1.2)
hold on
plot(tmin, h(1,:), '--', 'LineWidth',1.2)
plot(tmin, y(2,:), '-', 'LineWidth',1.2)
plot(tmin, h(2,:), '--', 'LineWidth',1.2)
xlabel('Time [min]')
ylabel('Liquid level [cm]')
title('Measured and True Tank Levels')
legend('Measured h_1','True h_1', ...
       'Measured h_2','True h_2','Location','best')
grid on

%% Save figure
plotdir = fullfile(here,'plots');
if ~exist(plotdir,'dir')
    mkdir(plotdir)
end

try, theme(gcf,'light'), catch, end
exportgraphics(gcf, ...
    fullfile(plotdir,'sim_2_3_sde.png'), ...
    'Resolution',150)
exportgraphics(gcf, ...
    fullfile(plotdir,'sim_2_3_sde.pdf'), ...
    'ContentType','vector')