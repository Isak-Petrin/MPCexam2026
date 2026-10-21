%% Modified Four Tank System - parameters and operating point
% Script, not a function: run it to define p, us, ds, xs, ys, zs and the
% noise parameters in the calling workspace.
% Physical constants from Week 2 slide 12, gamma from the exam text.

a1 = 1.2272;   % [cm^2] pipe cross sectional areas
a2 = 1.2272;
a3 = 1.2272;
a4 = 1.2272;

A1 = 380.1327; % [cm^2] tank cross sectional areas
A2 = 380.1327;
A3 = 380.1327;
A4 = 380.1327;

gamma1 = 0.60; % [-] valve constants (exam values)
gamma2 = 0.70;

g   = 981;     % [cm/s^2]
rho = 1.00;    % [g/cm^3]

p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];

ap = [a1; a2; a3; a4];
At = [A1; A2; A3; A4];

%% Operating point  (CHOICE)
us = [300; 300];   % [cm^3/s] nominal pump flows F1, F2  (Week 2 slide 35)
ds = [ 100;  80];   % [cm^3/s] nominal disturbance flows F3, F4

% Closed-form steady state (no toolbox needed). At steady state the
% outflow of each tank equals its total inflow and qout_i = a_i*sqrt(2*g*h_i),
% so h_i = (qout_i/a_i)^2/(2*g) and m_i = rho*A_i*h_i.
qs = zeros(4,1);
qs(3) = (1-gamma2)*us(2) + ds(1);      % tank 3: (1-gamma2)*F2 + F3
qs(4) = (1-gamma1)*us(1) + ds(2);      % tank 4: (1-gamma1)*F1 + F4
qs(1) = gamma1*us(1) + qs(3);          % tank 1: gamma1*F1 + outflow of tank 3
qs(2) = gamma2*us(2) + qs(4);          % tank 2: gamma2*F2 + outflow of tank 4
hs = (qs./ap).^2/(2*g);               % [cm]
xs = rho*At.*hs;                       % [g]

assert(norm(QuadrupleTankProcess(0,xs,us,ds,p)) < 1e-9, 'Steady state residual is not zero')
ys  = QuadrupleTankSensor(xs,p);
zs  = QuadrupleTankOutput(xs,p);

%% Noise parameters  (CHOICE)
Rdd   = diag([10 10].^2);    % [cm^3/s]^2   variance of piecewise constant d_k   (2.2)
Rvv   = diag([0.5 0.5].^2);  % [cm]^2       measurement noise variance           (2.2, 2.3)
sigma = [20; 20];            % [cm^3/s^0.5] diffusion of F3, F4                  (2.3)
seed  = 100;
