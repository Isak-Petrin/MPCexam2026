%% Exercise 2 - Stochastic 4-Tank System
% Continuous stochastic model:
%   dx(t) = A*x(t)dt + B( u(t)dt + sigma*dw(t) )
%         = (A*x(t) + B*u(t))dt + B*sigma*dw(t)
%
% sigma is 10% of the steady-state input u_s.
%
% Tasks:
% 1) Discretize the continuous stochastic system.
% 2) Simulate the discrete stochastic system.
% 3) Plot time series and histograms of the outputs.

clc
clear
close all

%% Parameters  % Slide 12

a1 = 1.2272;   % [cm^2]
a2 = 1.2272;   % [cm^2]
a3 = 1.2272;   % [cm^2]
a4 = 1.2272;   % [cm^2]

A1 = 380.1327; % [cm^2]
A2 = 380.1327; % [cm^2]
A3 = 380.1327; % [cm^2]
A4 = 380.1327; % [cm^2]

gamma1 = 0.45;
gamma2 = 0.40;

g   = 981;     % [cm/s^2]
rho = 1.00;    % [g/cm^3]

p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];


%% Steady state from Slide 35

us = [250; 325];   % [cm^3/s]

xs0 = 5000*ones(4,1);
xs = fsolve(@(x) QuadrupleTankProcess(0,x,us,p), xs0);

ys = QuadrupleTankSensor(xs,p);   % all four tank levels
zs = QuadrupleTankOutput(xs,p);   % controlled outputs [h1; h2]


%% Continuous-time linearization at steady state  % Slides 36-42

ap = [a1; a2; a3; a4];
At = [A1; A2; A3; A4];
gam = [gamma1; gamma2];

hs = ys;
Tconst = (At./ap).*sqrt(2*hs/g);

A = [-1/Tconst(1)   0               1/Tconst(3)    0;
      0              -1/Tconst(2)   0               1/Tconst(4);
      0               0             -1/Tconst(3)    0;
      0               0              0              -1/Tconst(4)];

B = [rho*gam(1)         0;
     0                  rho*gam(2);
     0                  rho*(1-gam(2));
     rho*(1-gam(1))     0];

% Controlled outputs z = [h1; h2]
C = [1/(rho*A1)  0             0  0;
     0            1/(rho*A2)   0  0];

D = zeros(2,2);


%% Continuous stochastic model  % Exercise 2 / Slide 70
%
% In deviation variables:
%
%   dX(t) = (A*X(t) + B*U(t))dt + G*dw(t)
%
% Exercise 2 gives
%
%   G = B*sigma
%
% where sigma is 10% of the steady-state input.

sigma = diag(0.10*us);
G = B*sigma;

fprintf('--- Exercise 2: Continuous stochastic model ---\n')
disp('sigma = 10% of steady-state inputs =')
disp(sigma)
disp('G = B*sigma =')
disp(G)


%% Discretization of the stochastic system  % Slides 68-70

% Use the same sampling time as Exercise 1.
Ts = 4;   % [s]

nx = size(A,1);
nu = size(B,2);

% -------------------------------------------------------------
% Deterministic part of Slide 70:
%
%   Abar = exp(A*Ts)
%   Bbar = integral_0^Ts exp(A*tau)*B dtau
%
% computed using the augmented matrix exponential.
% -------------------------------------------------------------

PhiAB = expm([A B;
              zeros(nu,nx) zeros(nu,nu)] * Ts);

Abar = PhiAB(1:nx,1:nx);
Bbar = PhiAB(1:nx,nx+1:nx+nu);

% -------------------------------------------------------------
% Stochastic part of Slide 70:
%
%   Gbar = I
%
%   Qbar = integral_0^Ts
%          exp(A*tau)*G*G'*exp(A'*tau) dtau
%
% with
%
%   [Phi11 Phi12;
%      0    Phi22]
%
%   = exp([-A   G*G';
%           0    A'] * Ts)
%
%   Qbar = Phi22' * Phi12
% -------------------------------------------------------------

Gbar = eye(nx);

PhiQ = expm([-A             G*G';
              zeros(nx,nx)  A'] * Ts);

Phi12 = PhiQ(1:nx,nx+1:2*nx);
Phi22 = PhiQ(nx+1:2*nx,nx+1:2*nx);

Qbar = Phi22' * Phi12;

% Remove tiny numerical asymmetry from floating-point arithmetic.
Qbar = (Qbar + Qbar')/2;

fprintf('\n--- Exercise 2: Discrete stochastic model ---\n')
disp('Abar =')
disp(Abar)
disp('Bbar =')
disp(Bbar)
disp('Gbar =')
disp(Gbar)
disp('Qbar =')
disp(Qbar)


%% Simulation of the discrete stochastic system  % Slide 70
%
%   X(k+1) = Abar*X(k) + Bbar*U(k) + Gbar*w(k)
%
%   w(k) ~ N(0,Qbar)
%
% We simulate around the steady state. Therefore:
%
%   X = x - xs
%   U = u - us
%
% and, since the physical input is held at u = us,
%
%   U = 0.

rng(100);   % reproducible random realization

t0 = 0;
tfinal = 20*60;   % [s]
t = t0:Ts:tfinal;
tmin = t/60;
N = length(t);

% Draw discrete process noise w_k ~ N(0,Qbar)
Lq = chol(Qbar,'lower');
w = Lq*randn(nx,N-1);

X = zeros(nx,N);
Zdev = zeros(2,N);

% Start exactly at the steady state: X(0) = 0
X(:,1) = zeros(nx,1);

U = zeros(nu,1);   % u = us -> U = 0

for k = 1:N-1
    Zdev(:,k) = C*X(:,k) + D*U;

    X(:,k+1) = Abar*X(:,k) ...
             + Bbar*U ...
             + Gbar*w(:,k);
end

Zdev(:,N) = C*X(:,N) + D*U;

% Convert deviation variables back to physical variables.
x = X + xs;
z = Zdev + zs;


%% Time-series plots of the outputs

figure;
subplot(2,1,1)
plot(tmin,z(1,:),'LineWidth',1.2)
yline(zs(1),'--','Steady state')
xlabel('time [min]')
ylabel('h_1 [cm]')
title('Stochastic output h_1')
grid on

subplot(2,1,2)
plot(tmin,z(2,:),'LineWidth',1.2)
yline(zs(2),'--','Steady state')
xlabel('time [min]')
ylabel('h_2 [cm]')
title('Stochastic output h_2')
grid on


%% Histograms of the outputs

figure;
subplot(1,2,1)
histogram(z(1,:))
xline(zs(1),'--','Steady state')
xlabel('h_1 [cm]')
ylabel('count')
title('Histogram of h_1')
grid on

subplot(1,2,2)
histogram(z(2,:))
xline(zs(2),'--','Steady state')
xlabel('h_2 [cm]')
ylabel('count')
title('Histogram of h_2')
grid on


%% Simple sample statistics of the simulated outputs

z_mean = mean(z,2);
z_std = std(z,0,2);

fprintf('\n--- Sample statistics from the stochastic simulation ---\n')
fprintf('h1: mean = %.4f cm, std = %.4f cm\n', z_mean(1), z_std(1));
fprintf('h2: mean = %.4f cm, std = %.4f cm\n', z_mean(2), z_std(2));
