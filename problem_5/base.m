%% 3 Continuous-time linearization at steady state  % Slides 36-42

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


%% 4 Discrete-time linear model, Ts = 4 s  % Slides 51-54

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


%% 5 Linear step responses and comparison with nonlinear model

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


%% 6 Continuous- and discrete-time transfer functions  % Slides 43-49 and 55-58

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