%% Problem 4.4 - transfer function identification (curve fitting of T3, T4)
% Identifies the 2x2 transfer function model
%
%            [ K11/(T1 s+1)               K12/((T1 s+1)(T3 s+1)) ]
%     G(s) = [                                                    ]
%            [ K21/((T2 s+1)(T4 s+1))     K22/(T2 s+1)            ]
%
% from the deterministic normalized step responses S_ij(t) = (h_i - h_i,s)/dF_j
% (same simulation as Problem 4.1).
%
%   1) Gains  K_ij = S_ij(t_end)                         (end value)
%   2) T1, T2      = time where S11, S22 reach 63.2% of K11, K22
%   3) T3, T4      = lsqcurvefit of the second-order step responses
%
%      S12_model(t;T3) = K12*(1 - (T1*exp(-t/T1) - T3*exp(-t/T3))/(T1 - T3))
%      S21_model(t;T4) = K21*(1 - (T2*exp(-t/T2) - T4*exp(-t/T4))/(T2 - T4))
%
%      T3_hat = argmin_{T3>0} sum_k (S12(t_k) - S12_model(t_k;T3))^2
%      T4_hat = argmin_{T4>0} sum_k (S21(t_k) - S21_model(t_k;T4))^2
%
% The fit is done for all three step sizes (+10/+25/+50 %); the +10 % step
% is used as the identified model (closest to the operating point).
% For comparison, the time constants of the linearized model,
% T_i = (A_i/a_i)*sqrt(2*h_i,s/g), are printed as well.
%
% Requires the Optimization Toolbox (lsqcurvefit).
% Results are saved to data/tf_identification.mat for Problems 4.5/4.6.

clc
clear
close all

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','problem_2','problem_2_1'));

QuadrupleTankParameters   % p, us, ds, xs, rho, At, ap, g, ...

%% Settings (same grid as Problem 4.1)
stepPct = [0.10 0.25 0.50];
Ts      = 10;                 % [s]
tfinal  = 20*60;              % [s]
t       = 0:Ts:tfinal;
tmin    = t/60;
nu      = 2;
ns      = length(stepPct);
jId     = 1;                  % step used for the identified model (+10 %)

hs = xs./(rho*At);            % steady-state levels [cm]

%% Deterministic step responses (as in sim_4_1_steps_deterministic.m)
S    = cell(nu,ns);           % S{i,j}: 4 x N normalized response, step in F_i
opts = odeset('RelTol',1e-8,'AbsTol',1e-6);
for i = 1:nu
    for j = 1:ns
        ustep    = us;
        du       = stepPct(j)*us(i);
        ustep(i) = us(i) + du;
        [~,X]    = ode15s(@(tt,xx) QuadrupleTankProcess(tt,xx,ustep,ds,p), t, xs, opts);
        S{i,j}   = (X'./(rho*At) - hs)/du;
    end
end

%% Model step responses
% First order: K*(1 - exp(-t/T))
% Second order (two first-order tanks in series, Ta ~= Tb):
%     K*(1 - (Ta*exp(-t/Ta) - Tb*exp(-t/Tb))/(Ta - Tb))
% with the critically damped limit used when Ta = Tb.
S2nd = @(t,K,Ta,Tb) secondOrderStep(t,K,Ta,Tb);

lsqopts = optimoptions('lsqcurvefit','Display','off', ...
                       'FunctionTolerance',1e-14,'StepTolerance',1e-12, ...
                       'OptimalityTolerance',1e-14);
T0 = 30;                      % [s] initial guess for T3, T4
lb = 1e-3;                    % T > 0
ub = Inf;

%% Identification for every step size
K     = zeros(2,2,ns);
T     = zeros(4,ns);
rmse  = zeros(2,ns);
fit12 = cell(1,ns);
fit21 = cell(1,ns);

for j = 1:ns
    S11 = S{1,j}(1,:);   S21 = S{1,j}(2,:);    % step in F1 -> h1, h2
    S12 = S{2,j}(1,:);   S22 = S{2,j}(2,:);    % step in F2 -> h1, h2

    % 1) gains = end values
    K(:,:,j) = [S11(end) S12(end);
                S21(end) S22(end)];

    % 2) T1, T2 from the 63.2 % rule (linear interpolation between samples)
    T1 = time63(t,S11);
    T2 = time63(t,S22);

    % 3) T3, T4 by least squares curve fitting
    K12 = K(1,2,j);
    K21 = K(2,1,j);
    [T3,res3] = lsqcurvefit(@(x,tt) S2nd(tt,K12,T1,x), T0, t, S12, lb, ub, lsqopts);
    [T4,res4] = lsqcurvefit(@(x,tt) S2nd(tt,K21,T2,x), T0, t, S21, lb, ub, lsqopts);

    T(:,j)    = [T1; T2; T3; T4];
    rmse(:,j) = sqrt([res3; res4]/numel(t));
    fit12{j}  = S2nd(t,K12,T1,T3);
    fit21{j}  = S2nd(t,K21,T2,T4);
end

%% Linearized model time constants for comparison
Tlin = (At./ap).*sqrt(2*hs/g);

%% Print results
fprintf('\nTransfer function identification (deterministic step responses)\n');
fprintf('%-6s %8s %8s %8s %8s | %7s %7s %7s %7s | %9s %9s\n', ...
        'step','K11','K12','K21','K22','T1','T2','T3','T4','RMSE S12','RMSE S21');
for j = 1:ns
    Kj = K(:,:,j);
    fprintf('+%3.0f%%  %8.4f %8.4f %8.4f %8.4f | %7.1f %7.1f %7.1f %7.1f | %9.2e %9.2e\n', ...
            100*stepPct(j),Kj(1,1),Kj(1,2),Kj(2,1),Kj(2,2),T(:,j),rmse(:,j));
end
fprintf('%-6s %8s %8s %8s %8s | %7.1f %7.1f %7.1f %7.1f |\n', ...
        'linear','','','','',Tlin);

%% Identified model (+10 % step)
Kid = K(:,:,jId);
Tid = T(:,jId);
s   = tf('s');
G   = [Kid(1,1)/(Tid(1)*s+1),                 Kid(1,2)/((Tid(1)*s+1)*(Tid(3)*s+1));
       Kid(2,1)/((Tid(2)*s+1)*(Tid(4)*s+1)),  Kid(2,2)/(Tid(2)*s+1)];
G.InputName  = {'F1','F2'};
G.OutputName = {'h1','h2'};
fprintf('\nIdentified model (+%.0f%% step):\n',100*stepPct(jId));
fprintf('  K  = [%.4f %.4f; %.4f %.4f]  [cm/(cm^3/s)]\n',Kid(1,1),Kid(1,2),Kid(2,1),Kid(2,2));
fprintf('  T1 = %.1f s, T2 = %.1f s, T3 = %.1f s, T4 = %.1f s\n',Tid);

%% Figure: data vs fitted model and residuals (+10 % step)
cols = [0.00 0.45 0.74; 0.85 0.33 0.10];
fig1 = figure('Name','4.4 curve fit T3 and T4','Position',[100 100 950 600]);

subplot(2,2,1)
plot(tmin,S{2,jId}(1,:),'o','Color',cols(1,:),'MarkerSize',4); hold on
plot(tmin,fit12{jId},'-','Color',cols(2,:),'LineWidth',1.4); hold off
grid on; box on; xlim(tmin([1 end]))
ylabel('\Deltah_1/\DeltaF_2  [cm/(cm^3/s)]')
title(sprintf('S_{12}: \\hat{T}_3 = %.1f s',Tid(3)))
legend('nonlinear simulation','fitted model','Location','southeast')

subplot(2,2,2)
plot(tmin,S{1,jId}(2,:),'o','Color',cols(1,:),'MarkerSize',4); hold on
plot(tmin,fit21{jId},'-','Color',cols(2,:),'LineWidth',1.4); hold off
grid on; box on; xlim(tmin([1 end]))
ylabel('\Deltah_2/\DeltaF_1  [cm/(cm^3/s)]')
title(sprintf('S_{21}: \\hat{T}_4 = %.1f s',Tid(4)))
legend('nonlinear simulation','fitted model','Location','southeast')

subplot(2,2,3)
plot(tmin,S{2,jId}(1,:) - fit12{jId},'k-','LineWidth',1.0)
grid on; box on; xlim(tmin([1 end]))
xlabel('time [min]'); ylabel('residual')

subplot(2,2,4)
plot(tmin,S{1,jId}(2,:) - fit21{jId},'k-','LineWidth',1.0)
grid on; box on; xlim(tmin([1 end]))
xlabel('time [min]'); ylabel('residual')

%% Figure: identified model vs simulation, all four channels (+10 % step)
fig2 = figure('Name','4.4 identified model vs simulation','Position',[150 100 950 600]);
Ystep = step(G,t);                        % N x 2 x 2  (time, output, input)
for i = 1:nu
    for r = 1:2
        subplot(2,2,(r-1)*nu+i)
        plot(tmin,S{i,jId}(r,:),'o','Color',cols(1,:),'MarkerSize',4); hold on
        plot(tmin,Ystep(:,r,i),'-','Color',cols(2,:),'LineWidth',1.4); hold off
        grid on; box on; xlim(tmin([1 end]))
        ylabel(sprintf('\\Deltah_%d/\\DeltaF_%d',r,i))
        if r == 1, title(sprintf('Step in F_%d',i)), end
        if r == 2, xlabel('time [min]'), end
        if r == 1 && i == nu
            legend('nonlinear simulation','identified G(s)','Location','southeast')
        end
    end
end

%% Save results
datadir = fullfile(here,'data');
if ~exist(datadir,'dir'), mkdir(datadir), end
save(fullfile(datadir,'tf_identification.mat'), ...
     'stepPct','K','T','rmse','Tlin','Kid','Tid','G','t');

plotdir = fullfile(here,'plots');
if ~exist(plotdir,'dir'), mkdir(plotdir), end
figs  = [fig1 fig2];
names = {'sim_4_4_curvefit_T3_T4','sim_4_4_identified_model'};
for k = 1:numel(figs)
    try, theme(figs(k),'light'), catch, end
    exportgraphics(figs(k),fullfile(plotdir,[names{k} '.png']),'Resolution',150)
    exportgraphics(figs(k),fullfile(plotdir,[names{k} '.pdf']),'ContentType','vector')
end

%% Local functions
function y = secondOrderStep(t,K,Ta,Tb)
% Step response of K/((Ta s+1)(Tb s+1))
if abs(Ta - Tb) < 1e-8*max(Ta,1)
    y = K*(1 - (1 + t/Ta).*exp(-t/Ta));
else
    y = K*(1 - (Ta*exp(-t/Ta) - Tb*exp(-t/Tb))/(Ta - Tb));
end
end

function T = time63(t,s)
% Time at which s first reaches 63.2 % of its final value
target = 0.632*s(end);
k = find(s >= target, 1, 'first');
T = t(k-1) + (target - s(k-1))*(t(k) - t(k-1))/(s(k) - s(k-1));
end
