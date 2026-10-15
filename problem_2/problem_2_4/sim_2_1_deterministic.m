%% Problem 2.4 - open-loop simulation of the deterministic model (2.1)

clc
clear
close all

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','problem_2_1'));

QuadrupleTankParameters   % p, us, ds, xs, ys, zs, At, rho, ...

%% Time grid
Ts     = 10;        % [s] sampling time of the simulation grid
tfinal = 20*60;     % [s]
t      = 0:Ts:tfinal;
tmin   = t/60;
N      = length(t);

nx = 4;   % states
ny = 2;   % measured levels h1, h2
nz = 2;   % controlled levels h1, h2

%% Piecewise constant input sequences (Week 2 slides, one value per minute)
F1 = [330 235 295 303 242 292 300 243 325 324 ...
      330 357 278 274 202 301 320 334 276 308];

F2 = [285 301 295 313 314 316 289 278 302 293 ...
      304 296 307 277 319 305 300 306 307 297];

U = [F1; F2];
u = [repelem(U,1,60/Ts) U(:,end)];   % hold each value for 60 s -> 2 x N

% Disturbances: nominal values with a step in F3 after 10 min  (CHOICE)
d = repmat(ds,1,N);
d(1, t >= 10*60) = ds(1) + 50;

%% Simulation: ode15s over each sampling interval
x = zeros(nx,N);
y = zeros(ny,N);
z = zeros(nz,N);
x(:,1) = xs;

for k = 1:N-1
    y(:,k) = QuadrupleTankSensor(x(:,k),p);
    z(:,k) = QuadrupleTankOutput(x(:,k),p);

    [~,Xk] = ode15s(@QuadrupleTankProcess, ...
                    [t(k) t(k+1)], ...
                    x(:,k), ...
                    [], ...
                    u(:,k), d(:,k), p);

    x(:,k+1) = Xk(end,:)';
end

y(:,N) = QuadrupleTankSensor(x(:,N),p);
z(:,N) = QuadrupleTankOutput(x(:,N),p);

h = x./(rho*At);   % all four levels [cm], for plotting only

%% Plots: inputs and disturbances
figure('Name','2.1 inputs');
subplot(2,2,1)
stairs(tmin,u(1,:),'LineWidth',1.5)
xlabel('time [min]'); ylabel('F_1 [cm^3/s]'); title('Pump flow F_1'); grid on
subplot(2,2,2)
stairs(tmin,u(2,:),'LineWidth',1.5)
xlabel('time [min]'); ylabel('F_2 [cm^3/s]'); title('Pump flow F_2'); grid on
subplot(2,2,3)
stairs(tmin,d(1,:),'LineWidth',1.5)
xlabel('time [min]'); ylabel('F_3 [cm^3/s]'); title('Disturbance F_3'); grid on
subplot(2,2,4)
stairs(tmin,d(2,:),'LineWidth',1.5)
xlabel('time [min]'); ylabel('F_4 [cm^3/s]'); title('Disturbance F_4'); grid on

%% Plots: outputs
figure('Name','2.1 outputs');
subplot(2,2,1)
plot(tmin,h(3,:),'LineWidth',1.2)
xlabel('time [min]'); ylabel('h_3 [cm]'); title('Tank 3'); grid on
subplot(2,2,2)
plot(tmin,h(4,:),'LineWidth',1.2)
xlabel('time [min]'); ylabel('h_4 [cm]'); title('Tank 4'); grid on
subplot(2,2,3)
plot(tmin,y(1,:),'LineWidth',1.2)
xlabel('time [min]'); ylabel('h_1 [cm]'); title('Tank 1 (y_1 = z_1)'); grid on
subplot(2,2,4)
plot(tmin,y(2,:),'LineWidth',1.2)
xlabel('time [min]'); ylabel('h_2 [cm]'); title('Tank 2 (y_2 = z_2)'); grid on
