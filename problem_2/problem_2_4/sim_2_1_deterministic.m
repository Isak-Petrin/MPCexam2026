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
% d(1, t >= 10*60) = ds(1) + 50;

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

%% Plots: inputs and outputs
figure('Name','2.1 deterministic');

subplot(3,1,1)
stairs(tmin,u','LineWidth',1.2)
ylabel('[cm^3/s]'); title('Pump flows u = [F_1; F_2]')
legend('F_1','F_2','Location','best'); grid on

subplot(3,1,2)
stairs(tmin,d','LineWidth',1.2)
ylabel('[cm^3/s]'); title('Disturbance flows d = [F_3; F_4]')
legend('F_3','F_4','Location','best'); grid on

subplot(3,1,3)
plot(tmin,y','LineWidth',1.2)
xlabel('time [min]'); ylabel('[cm]'); title('Measured levels y = [h_1; h_2]')
legend('h_1','h_2','Location','best'); grid on

%% Save figure to problem_2_4/plots
plotdir = fullfile(here,'plots');
if ~exist(plotdir,'dir')
    mkdir(plotdir)
end
try, theme(gcf,'light'), catch, end   % report figures on white (R2025a+)
exportgraphics(gcf,fullfile(plotdir,'sim_2_1_deterministic.png'),'Resolution',150)
exportgraphics(gcf,fullfile(plotdir,'sim_2_1_deterministic.pdf'),'ContentType','vector')
