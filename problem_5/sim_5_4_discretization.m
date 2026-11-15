%% Problem 5.4 - discrete-time linear state space models
% Linearizes the modified four tank system at the operating point from
% QuadrupleTankParameters (Section 5.1) and computes the ZOH
% discretization for the three models from Problem 2:
%
%   deterministic           X_{k+1} = Ad X_k + Bd U_k + Ed D_k
%   piecewise constant d_k  as above with D_k ~ N(0,Rdd), v_k ~ N(0,Rvv)
%   SDE                     as above with w_k ~ N(0,Qd),  v_k ~ N(0,Rvv)
%
% The result is verified by numerical quadrature, against exp(-Ts/T_i) and by a
% continuous vs. discrete step response. Everything is saved to
% data/linear_models.mat for Problems 5.5, 6 and 8-11.

clc
clear
close all

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','problem_2','problem_2_1'));

QuadrupleTankParameters   % p, us, ds, xs, ys, zs, rho, At, ap, sigma, Rdd, Rvv

%% Settings
Ts = 10;                  % [s] sampling time (T3 ~ 49 s -> ~5 samples per time constant)

%% Continuous-time linearization (Section 5.1)
hs = xs./(rho*At);
T  = (At./ap).*sqrt(2*hs/g);              % time constants (5.13)

A = [-1/T(1)  0       1/T(3)  0;
      0      -1/T(2)  0       1/T(4);
      0       0      -1/T(3)  0;
      0       0       0      -1/T(4)];

B = [rho*gamma1        0;
     0                 rho*gamma2;
     0                 rho*(1-gamma2);
     rho*(1-gamma1)    0];

E = [0 0; 0 0; rho 0; 0 rho];

% Sensor and output matrices: the functions are linear, so a finite
% difference gives the exact Jacobian and stays consistent with the code.
nx = 4;
ny = length(ys);
nz = length(zs);
C  = zeros(ny,nx);
Cz = zeros(nz,nx);
for j = 1:nx
    dx = zeros(nx,1); dx(j) = 1;
    C(:,j)  = QuadrupleTankSensor(xs+dx,p) - ys;
    Cz(:,j) = QuadrupleTankOutput(xs+dx,p) - zs;
end

% Diffusion matrix of the SDE model (5.19)
G = E*diag(sigma);

%% Check analytical Jacobians by finite differences (5.21)
fx = @(x) QuadrupleTankProcess(0,x,us,ds,p);
fu = @(u) QuadrupleTankProcess(0,xs,u,ds,p);
fd = @(d) QuadrupleTankProcess(0,xs,us,d,p);
Afd = zeros(nx,nx); Bfd = zeros(nx,2); Efd = zeros(nx,2);
for j = 1:nx
    del = sqrt(eps)*(1+abs(xs(j))); e = zeros(nx,1); e(j) = del;
    Afd(:,j) = (fx(xs+e)-fx(xs))/del;
end
for j = 1:2
    del = sqrt(eps)*(1+abs(us(j))); e = zeros(2,1); e(j) = del;
    Bfd(:,j) = (fu(us+e)-fu(us))/del;
    del = sqrt(eps)*(1+abs(ds(j))); e = zeros(2,1); e(j) = del;
    Efd(:,j) = (fd(ds+e)-fd(ds))/del;
end
fprintf('Jacobian check: max|A-Afd| = %.1e, max|B-Bfd| = %.1e, max|E-Efd| = %.1e\n', ...
        max(abs(A-Afd),[],'all'), max(abs(B-Bfd),[],'all'), max(abs(E-Efd),[],'all'));

%% Discretization (Problem 5.4)
[Ad,Bd,Ed] = DiscretizeZOH(A,B,E,Ts);     % deterministic part
[AdVL,Qd]  = DiscretizeSDE(A,G,Ts);       % process noise of the SDE model
Qdd        = Ed*Rdd*Ed';                  % state noise from piecewise constant d_k (2.2)

fprintf('\nDiscrete-time model, Ts = %g s\n',Ts);
disp('Ad =');  disp(Ad)
disp('Bd =');  disp(Bd)
disp('Ed =');  disp(Ed)
disp('Qd (SDE) =');           disp(Qd)
disp('Ed*Rdd*Ed'' (2.2) =');  disp(Qdd)

%% Verification
% 1) against direct numerical integration of (5.x) (no toolbox needed)
Gam = integral(@(s) expm(A*s), 0, Ts, 'ArrayValued', true);   % int_0^Ts e^{As} ds
QdI = integral(@(s) expm(A*s)*(G*G')*expm(A'*s), 0, Ts, 'ArrayValued', true);
fprintf('Quadrature check: max|[Bd Ed]-Gam*[B E]| = %.1e, max|Qd-QdI| = %.1e, max|Ad-Ad_VanLoan| = %.1e\n', ...
        max(abs([Bd Ed]-Gam*[B E]),[],'all'), max(abs(Qd-QdI),[],'all'), max(abs(Ad-AdVL),[],'all'));

% 2) eigenvalues of Ad equal exp(-Ts/T_i)
lam = sort(eig(Ad));
fprintf('\n%6s %12s %12s\n','tank','eig(Ad)','exp(-Ts/Ti)');
lamT = exp(-Ts./T);
[~,idx] = sort(lamT);
for i = 1:nx
    fprintf('%6d %12.4f %12.4f\n', idx(i), lam(i), lamT(idx(i)));
end

% 3) steady-state gain of the discrete model equals G(0)
Kd = Cz*((eye(nx)-Ad)\Bd);
Kc = -Cz*(A\B);
fprintf('\nGain check: max|Cz(I-Ad)^{-1}Bd - (-Cz A^{-1} B)| = %.1e\n', max(abs(Kd-Kc),[],'all'));
disp('K ='); disp(Kd)

%% Continuous vs. discrete step response (+10% step in F1 and in F2)
tfinal = 20*60;
tc     = linspace(0,tfinal,2001);        % fine grid for the continuous model
tk     = 0:Ts:tfinal;                    % sampling instants
Nk     = length(tk);

Zc = cell(1,2); Zk = cell(1,2);
for i = 1:2
    U = zeros(2,1); U(i) = 0.10*us(i);

    % continuous: integrate Xdot = A X + B U with ode45 (constant input)
    [~,Xc] = ode45(@(tt,xx) A*xx + B*U, tc, zeros(nx,1), odeset('RelTol',1e-10,'AbsTol',1e-10));
    Zc{i}  = Cz*Xc';

    % discrete: the recursion X_{k+1} = Ad X_k + Bd U_k
    X = zeros(nx,1); Zk{i} = zeros(nz,Nk);
    for k = 1:Nk
        Zk{i}(:,k) = Cz*X;
        X = Ad*X + Bd*U;
    end
end

err = max(cellfun(@(zc,zk) max(abs(interp1(tc,zc',tk)'-zk),[],'all'), Zc, Zk));
fprintf('Max |continuous - discrete| at the sampling instants: %.1e cm\n', err);

fig1 = figure('Name','5.4 continuous vs discrete','Position',[100 100 900 600]);
outname = {'h_1','h_2'};
for i = 1:2
    for r = 1:nz
        subplot(nz,2,(r-1)*2+i)
        plot(tc/60, Zc{i}(r,:),'LineWidth',1.3); hold on
        plot(tk/60, Zk{i}(r,:),'o','MarkerSize',3.5); hold off
        grid on; box on; xlim([0 tfinal/60])
        ylabel(['\Delta' outname{r} ' [cm]'])
        if r == 1, title(sprintf('+10%% step in F_%d',i)), end
        if r == nz, xlabel('time [min]'), end
        if r == 1 && i == 2
            legend('continuous','discrete, T_s = 10 s','Location','southeast')
        end
    end
end

%% Save
datadir = fullfile(here,'data');
if ~exist(datadir,'dir'), mkdir(datadir), end
save(fullfile(datadir,'linear_models.mat'), ...
     'Ts','A','B','E','C','Cz','G','T','Ad','Bd','Ed','Qd','Qdd','Rdd','Rvv', ...
     'us','ds','xs','ys','zs','hs');

plotdir = fullfile(here,'plots');
if ~exist(plotdir,'dir'), mkdir(plotdir), end
try, theme(fig1,'light'), catch, end
exportgraphics(fig1,fullfile(plotdir,'sim_5_4_continuous_vs_discrete.png'),'Resolution',150)
exportgraphics(fig1,fullfile(plotdir,'sim_5_4_continuous_vs_discrete.pdf'),'ContentType','vector')
