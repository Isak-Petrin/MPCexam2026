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

nx = 4;
ny = 2;
nz = 2;
nw = 2;             % Wiener processes: one for F3, one for F4

%% Piecewise constant input sequences (same as sim_2_1_deterministic)
F1 = [330 235 295 303 242 292 300 243 325 324 ...
      330 357 278 274 202 301 320 334 276 308];

F2 = [285 301 295 313 314 316 289 278 302 293 ...
      304 296 307 277 319 305 300 306 307 297];

U = [F1; F2];
u = [repelem(U,1,60/Ts) U(:,end)];

% Mean disturbance with a step in F3 after 10 min (CHOICE).
% The fluctuations around the mean come from sigma*dw in the SDE.
d = repmat(ds,1,N);
d(1, t >= 10*60) = ds(1) + 50;

%% Wiener process for the whole horizon (nw channels, one realization)
[W,Tw,~] = StdWienerProcess(tfinal,(N-1)*nsub,nw,1,seed);   % W: nw x ((N-1)*nsub+1)

%% Measurement noise
Lr = chol(Rvv,'lower');
v  = Lr*randn(ny,N);

%% Simulation: Euler-Maruyama over each sampling interval
psde = [p; sigma];   % p(13:14) = sigma3, sigma4 for QuadrupleTankDiffusion

x = zeros(nx,N);
y = zeros(ny,N);
z = zeros(nz,N);
x(:,1) = xs;

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

%% Plots: inputs and mean disturbances
figure('Name','2.3 inputs');
subplot(2,2,1)
stairs(tmin,u(1,:),'LineWidth',1.5)
xlabel('time [min]'); ylabel('F_1 [cm^3/s]'); title('Pump flow F_1'); grid on
subplot(2,2,2)
stairs(tmin,u(2,:),'LineWidth',1.5)
xlabel('time [min]'); ylabel('F_2 [cm^3/s]'); title('Pump flow F_2'); grid on
subplot(2,2,3)
stairs(tmin,d(1,:),'LineWidth',1.5)
xlabel('time [min]'); ylabel('F_3 [cm^3/s]'); title('Mean disturbance F_3'); grid on
subplot(2,2,4)
stairs(tmin,d(2,:),'LineWidth',1.5)
xlabel('time [min]'); ylabel('F_4 [cm^3/s]'); title('Mean disturbance F_4'); grid on

%% Plots: outputs, noisy measurement vs true level
figure('Name','2.3 outputs');
subplot(2,2,1)
plot(tmin,h(3,:),'LineWidth',1.2)
xlabel('time [min]'); ylabel('h_3 [cm]'); title('Tank 3 (unmeasured)'); grid on
subplot(2,2,2)
plot(tmin,h(4,:),'LineWidth',1.2)
xlabel('time [min]'); ylabel('h_4 [cm]'); title('Tank 4 (unmeasured)'); grid on
subplot(2,2,3)
plot(tmin,y(1,:),'-',tmin,h(1,:),'--','LineWidth',1.2)
xlabel('time [min]'); ylabel('h_1 [cm]'); title('Tank 1')
legend('Noisy measurement y_1','True level','Location','best'); grid on
subplot(2,2,4)
plot(tmin,y(2,:),'-',tmin,h(2,:),'--','LineWidth',1.2)
xlabel('time [min]'); ylabel('h_2 [cm]'); title('Tank 2')
legend('Noisy measurement y_2','True level','Location','best'); grid on

%% Save figure to problem_2_4/plots
plotdir = fullfile(here,'plots');
if ~exist(plotdir,'dir')
    mkdir(plotdir)
end
try, theme(gcf,'light'), catch, end   % report figures on white (R2025a+)
exportgraphics(gcf,fullfile(plotdir,'sim_2_3_sde.png'),'Resolution',150)
exportgraphics(gcf,fullfile(plotdir,'sim_2_3_sde.pdf'),'ContentType','vector')
