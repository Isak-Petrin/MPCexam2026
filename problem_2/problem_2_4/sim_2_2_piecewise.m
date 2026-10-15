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

% Nominal (mean) disturbance with a step in F3 after 10 min  (CHOICE)
dnom = repmat(ds,1,N);
dnom(1, t >= 10*60) = ds(1) + 50;

%% Noise realizations
rng(seed);

Lq = chol(Rdd,'lower');
w  = Lq*randn(2,N);          % disturbance noise, one draw per sampling interval

Lr = chol(Rvv,'lower');
v  = Lr*randn(ny,N);         % measurement noise

d = dnom + w;                % realized piecewise constant disturbance d_k

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
z(:,N) = QuadrupleTankOutput(x(:,N),p);

h = x./(rho*At);   % true levels [cm]

%% Plots: inputs and realized disturbances
figure('Name','2.2 inputs');
subplot(2,2,1)
stairs(tmin,u(1,:),'LineWidth',1.5)
xlabel('time [min]'); ylabel('F_1 [cm^3/s]'); title('Pump flow F_1'); grid on
subplot(2,2,2)
stairs(tmin,u(2,:),'LineWidth',1.5)
xlabel('time [min]'); ylabel('F_2 [cm^3/s]'); title('Pump flow F_2'); grid on
subplot(2,2,3)
stairs(tmin,d(1,:),'LineWidth',1.2); hold on
stairs(tmin,dnom(1,:),'--','LineWidth',1.2)
xlabel('time [min]'); ylabel('F_3 [cm^3/s]'); title('Disturbance F_3')
legend('realized d_k','mean','Location','best'); grid on
subplot(2,2,4)
stairs(tmin,d(2,:),'LineWidth',1.2); hold on
stairs(tmin,dnom(2,:),'--','LineWidth',1.2)
xlabel('time [min]'); ylabel('F_4 [cm^3/s]'); title('Disturbance F_4')
legend('realized d_k','mean','Location','best'); grid on

%% Plots: outputs, noisy measurement vs true level
figure('Name','2.2 outputs');
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
