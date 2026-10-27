% Week 2 "Slide recreations.m" lines 177-310

clc
clear
close all

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','problem_2_1'));

QuadrupleTankParameters   % p, us, ds, xs, Rdd, Rvv, seed, ...

%% Time grid
Ts     = 10;        % [s]
tfinal = 20*60;     % [s]
t      = 0:Ts:tfinal;
tmin   = t/60;
N      = length(t);

nx = 4;
ny = 2;
nz = 2;

%% Piecewise constant input sequences (same as sim_2_1_deterministic)
F1 = [330 235 295 303 242 292 300 243 325 324 ...
      330 357 278 274 202 301 320 334 276 308];

F2 = [285 301 295 313 314 316 289 278 302 293 ...
      304 296 307 277 319 305 300 306 307 297];

U = [F1; F2];
u = [repelem(U,1,60/Ts) U(:,end)];

rng(seed);

% Nominal (mean) disturbance with a step in F3 after 10 min  (CHOICE)
dnom = repmat(ds,1,N);
%dnom(1, t >= 10*60) = ds(1) + 50; %STEP INCREASE

%% Disturbances: d_i = [d_i(t_1) ... d_i(t_N)]' ~ N(mu_i, R_d),  R_d in R^{N x N}
rng(seed);

ell = 2*60;                          % correlation length of R_d [s]  (CHOICE)
[T1,T2] = meshgrid(t,t);
K = eye(N);                          % N x N identity -> independent samples in time    % N x N temporal correlation (squared exponential)

d = zeros(2,N);
for i = 1:2
    mu_i = dnom(i,:)';               % mean vector mu in R^N
    Rd_i = Rdd(i,i)*K;               % covariance R_d in R^{N x N}
    % Factor R_d = L*L'. chol() fails here: the SE kernel is numerically singular.
    [V,D] = eig((Rd_i+Rd_i')/2);
    Ld = V*diag(sqrt(max(diag(D),0)));
    d(i,:) = (mu_i + Ld*randn(N,1))';    % one draw of the whole vector d_i
end
w = d - dnom;                        % deviation from mean (kept for reference)

%% Measurement noise
vscale = 20
Lr = chol(vscale * Rvv,'lower');
v  = Lr*randn(ny,N);           % realized piecewise constant disturbance d_k

%% Simulation: ode15s over each sampling interval, d_k held constant
x = zeros(nx,N);
y = zeros(ny,N);
z = zeros(nz,N);
x(:,1) = xs;

for k = 1:N-1
    y(:,k) = QuadrupleTankSensor(x(:,k),p) + v(:,k);   % y = g(x,p) + v
    z(:,k) = QuadrupleTankOutput(x(:,k),p);

    [~,Xk] = ode15s(@QuadrupleTankProcess, ...
                    [t(k) t(k+1)], ...
                    x(:,k), ...
                    [], ...
                    u(:,k), d(:,k), p);

    x(:,k+1) = Xk(end,:)';
end

y(:,N) = QuadrupleTankSensor(x(:,N),p) + v(:,N);
z(:,N) = QuadrupleTankOutput(x(:,N),p) + v(1:2,N);

h = x./(rho*At);   % true levels [cm]

%% Plots: inputs and outputs
figure('Name','2.2 piecewise constant stochastic');

subplot(3,1,1)
stairs(tmin,u','LineWidth',1.2)
ylabel('[cm^3/s]'); title('Pump flows u = [F_1; F_2]')
legend('F_1','F_2','Location','best'); grid on

subplot(3,1,2)
stairs(tmin,d','LineWidth',1.2)
hold on
stairs(tmin,dnom','k--','LineWidth',1.0)
ylabel('[cm^3/s]'); title('Disturbance flows d_k = mean + w_k')
legend('F_3','F_4','mean','Location','best'); grid on

subplot(3,1,3)
plot(tmin,y','LineWidth',1.2)
hold on
plot(tmin,h(1:2,:)','k--','LineWidth',1.0)
xlabel('time [min]'); ylabel('[cm]'); title('Measured levels y = g(x,p) + v')
legend('y_1','y_2','true h_1, h_2','Location','best'); grid on

%% Save figure to problem_2_4/plots
plotdir = fullfile(here,'plots');
if ~exist(plotdir,'dir')
    mkdir(plotdir)
end
try, theme(gcf,'light'), catch, end   % report figures on white (R2025a+)
exportgraphics(gcf,fullfile(plotdir,'sim_2_2_piecewise.png'),'Resolution',150)
exportgraphics(gcf,fullfile(plotdir,'sim_2_2_piecewise.pdf'),'ContentType','vector')
