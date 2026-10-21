%% Exercise 1 - 4-Tank System
% F1s = 250 cm^3/s, F2s = 325 cm^3/s
% 1) Compute steady state
% 2) Nonlinear responses to 5%, 10%, 25% step in F1
% 3) Linearize at steady state -> continuous-time (A,B,C,D)
% 4) Discretize with Ts = 4 s -> discrete-time (Ad,Bd,Cd,Dd)
% 5) Linear step responses and comparison with nonlinear model
% 6) Continuous/discrete transfer functions, gains, poles, zeros

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


%% 1) Steady state  % Slide 35

us = [250; 325]; % [cm^3/s]

xs0 = 5000*ones(4,1);
xs = fsolve(@(x) QuadrupleTankProcess(0,x,us,p), xs0);

ys = QuadrupleTankSensor(xs,p);
zs = QuadrupleTankOutput(xs,p);

disp('--- Exercise 1.1: Steady state ---')
disp('us = [F1s; F2s] =')
disp(us)
disp('xs = [m1s; m2s; m3s; m4s] [g] =')
disp(xs)
disp('ys = [h1s; h2s; h3s; h4s] [cm] =')
disp(ys)
disp('zs = [h1s; h2s] [cm] =')
disp(zs)


%% 2) Nonlinear step responses: +5%, +10%, +25% in F1  % Slides 18-24

stepPct = [0.05 0.10 0.25];
Ts_sim = 4;          % [s], use same grid as Exercise 1.4
tfinal = 20*60;     % [s]
t = 0:Ts_sim:tfinal;
tmin = t/60;
N = length(t);

nx = 4;
nz = 2;

xNL = cell(size(stepPct));
zNL = cell(size(stepPct));

for j = 1:length(stepPct)
    x = zeros(nx,N);
    z = zeros(nz,N);
    x(:,1) = xs;

    ustep = [us(1)*(1+stepPct(j)); us(2)];

    for k = 1:N-1
        z(:,k) = QuadrupleTankOutput(x(:,k),p);

        [~,Xk] = ode15s(@QuadrupleTankProcess, ...
                        [t(k) t(k+1)], ...
                        x(:,k), ...
                        [], ...
                        ustep, ...
                        p);

        x(:,k+1) = Xk(end,:)';
    end

    z(:,N) = QuadrupleTankOutput(x(:,N),p);

    xNL{j} = x;
    zNL{j} = z;
end

figure;
for j = 1:length(stepPct)
    subplot(3,1,j)
    plot(tmin,zNL{j}(1,:),'LineWidth',1.2)
    hold on
    plot(tmin,zNL{j}(2,:),'LineWidth',1.2)
    yline(zs(1),'--')
    yline(zs(2),'--')
    xlabel('time [min]')
    ylabel('level [cm]')
    title(sprintf('Nonlinear response: +%.0f%% step in F_1',100*stepPct(j)))
    legend('h_1','h_2','h_{1,s}','h_{2,s}','Location','best')
    grid on
end


%% 3) Continuous-time linearization at steady state  % Slides 36-42

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

% Exercise asks for the 2 controlled outputs z = [h1; h2]
C = [1/(rho*A1)  0             0  0;
     0            1/(rho*A2)   0  0];

D = zeros(2,2);

disp('--- Exercise 1.3: Continuous-time linear model ---')
disp('A ='); disp(A)
disp('B ='); disp(B)
disp('C ='); disp(C)
disp('D ='); disp(D)

sys_c = ss(A,B,C,D);


%% 4) Discrete-time linear model, Ts = 4 s  % Slides 51-54

Ts = 4; % [s]
nu = size(B,2);

M = expm([A B;
          zeros(nu,nx) zeros(nu,nu)] * Ts);

Ad = M(1:nx,1:nx);
Bd = M(1:nx,nx+1:nx+nu);
Cd = C;
Dd = D;

disp('--- Exercise 1.4: Discrete-time linear model, Ts = 4 s ---')
disp('Ad ='); disp(Ad)
disp('Bd ='); disp(Bd)
disp('Cd ='); disp(Cd)
disp('Dd ='); disp(Dd)

sys_d = ss(Ad,Bd,Cd,Dd,Ts);


%% 5) Linear step responses and comparison with nonlinear model

zLIN = cell(size(stepPct));

for j = 1:length(stepPct)
    X = zeros(nx,N);
    Zdev = zeros(nz,N);

    % Step in deviation input U = u-us
    Ustep = [stepPct(j)*us(1); 0];

    for k = 1:N-1
        Zdev(:,k) = C*X(:,k) + D*Ustep;
        X(:,k+1) = Ad*X(:,k) + Bd*Ustep;
    end
    Zdev(:,N) = C*X(:,N) + D*Ustep;

    zLIN{j} = Zdev + zs;
end

figure;
for j = 1:length(stepPct)
    subplot(3,2,2*j-1)
    plot(tmin,zNL{j}(1,:),'LineWidth',1.2)
    hold on
    plot(tmin,zLIN{j}(1,:),'--','LineWidth',1.2)
    xlabel('time [min]')
    ylabel('h_1 [cm]')
    title(sprintf('Tank 1: +%.0f%% F_1 step',100*stepPct(j)))
    legend('Nonlinear','Linear','Location','best')
    grid on

    subplot(3,2,2*j)
    plot(tmin,zNL{j}(2,:),'LineWidth',1.2)
    hold on
    plot(tmin,zLIN{j}(2,:),'--','LineWidth',1.2)
    xlabel('time [min]')
    ylabel('h_2 [cm]')
    title(sprintf('Tank 2: +%.0f%% F_1 step',100*stepPct(j)))
    legend('Nonlinear','Linear','Location','best')
    grid on
end

% Simple numerical comparison: maximum absolute output mismatch
for j = 1:length(stepPct)
    err = max(abs(zNL{j} - zLIN{j}),[],2);
    fprintf('Max |nonlinear-linear| for +%.0f%% F1 step: h1 = %.4f cm, h2 = %.4f cm\n', ...
            100*stepPct(j), err(1), err(2));
end


%% 6) Continuous- and discrete-time transfer functions  % Slides 43-49 and 55-58

Gc = tf(sys_c);
Gd = tf(sys_d);

disp('--- Exercise 1.6: Continuous-time transfer function matrix ---')
Gc

disp('--- Exercise 1.6: Discrete-time transfer function matrix ---')
Gd

% Steady-state gain matrices
Kc = dcgain(sys_c);
Kd = dcgain(sys_d);

disp('Continuous-time steady-state gain Kc = Gc(0) =')
disp(Kc)

disp('Discrete-time steady-state gain Kd = Gd(1) =')
disp(Kd)

% State-space poles
poles_c = pole(sys_c);
poles_d = pole(sys_d);

disp('Continuous-time poles =')
disp(poles_c)

disp('Discrete-time poles =')
disp(poles_d)

% Transmission zeros of the full 2x2 MIMO systems
zeros_c = tzero(sys_c);
zeros_d = tzero(sys_d);

disp('Continuous-time transmission zeros =')
disp(zeros_c)

disp('Discrete-time transmission zeros =')
disp(zeros_d)

% For completeness, print poles/zeros of each SISO channel as well.
for iy = 1:2
    for iu = 1:2
        Gc_ij = Gc(iy,iu);
        Gd_ij = Gd(iy,iu);

        fprintf('\nChannel y%d/u%d\n',iy,iu)
        fprintf('Continuous gain = %.8g\n',dcgain(Gc_ij))
        disp('Continuous poles ='); disp(pole(Gc_ij))
        disp('Continuous zeros ='); disp(zero(Gc_ij))

        fprintf('Discrete gain = %.8g\n',dcgain(Gd_ij))
        disp('Discrete poles ='); disp(pole(Gd_ij))
        disp('Discrete zeros ='); disp(zero(Gd_ij))
    end
end
