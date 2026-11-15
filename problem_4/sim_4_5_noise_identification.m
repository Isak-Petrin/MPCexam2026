%% Problem 4.5 - identification from noisy step responses (Monte Carlo)
% Repeats the transfer function identification on the stochastic step
% responses (SDE model + measurement noise, as in
% sim_4_2_steps_stochastic.m) for many noise realizations, to see how the
% identified gains and time constants depend on the noise level and on
% the step size.
%
% For each noise level (low/medium/high), step size (+10/+25/+50%) and
% noise realization, the model
%
%            [ K11/(T1 s+1)               K12/((T1 s+1)(T3 s+1)) ]
%     G(s) = [                                                    ]
%            [ K21/((T2 s+1)(T4 s+1))     K22/(T2 s+1)            ]
%
% is identified from the normalized measured step responses
%     S_ij(t_k) = (y_i(t_k) - y_i,s)/dF_j
% with a method that does not rely on single samples:
%   1) gains  K_ij = mean of S_ij(t_k) over the settled part, t >= tTail
%   2) T1, T2      = least-squares fit of the first-order step response
%                    to S11, S22 (gain fixed)
%   3) T3, T4      = least-squares fit of the second-order step response
%                    to S12, S21 (gain and T1, T2 fixed), as in sim_4_4
% The same method is applied to the deterministic step responses, and all
% estimates are compared with the linearized model at the operating point.
%
% Noise: realization 1 uses the same random numbers as
% sim_4_2_steps_stochastic.m (seed and seed+1). Within one realization the
% same Wiener path and measurement noise are used for all noise levels,
% step sizes and inputs.
%
% Runs in base MATLAB (fminbnd; no toolboxes).
%
% Output:
%   - printed tables: mean and standard deviation of every parameter, and
%     the relative RMS error with respect to the linearized model
%   - plots/sim_4_5_noise_identification.{png,pdf}: relative RMS error
%     against step size for each noise level
%   - data/noise_identification.mat

clc
clear
close all

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','problem_2','problem_2_1'));
addpath(fullfile(here,'..','problem_2','problem_2_4'));

QuadrupleTankParameters   % p, us, ds, xs, ys, rho, At, ap, g, gamma1, gamma2, seed, ...

%% Settings (simulation as in sim_4_2_steps_stochastic.m)
stepPct = [0.10 0.25 0.50];
Ts      = 10;            % [s] sampling time (measurements)
tfinal  = 20*60;         % [s]
t       = 0:Ts:tfinal;
N       = length(t);
nsub    = 10;            % Euler-Maruyama sub-steps per sample (dt = 1 s)
Nmc     = 500;           % number of noise realizations
tTail   = 10*60;         % [s] gains are averaged over t >= tTail

nu = 2;
ny = 2;
nw = 2;
ns = length(stepPct);

% Noise levels (same as sim_4_2_steps_stochastic.m)
noise(1).name  = 'low';
noise(1).sigma = [ 5;  5];            % [cm^3/s^0.5] diffusion on F3, F4
noise(1).Rvv   = diag([0.1 0.1].^2);  % [cm^2]

noise(2).name  = 'medium';
noise(2).sigma = [20; 20];
noise(2).Rvv   = diag([0.5 0.5].^2);

noise(3).name  = 'high';
noise(3).sigma = [50; 50];
noise(3).Rvv   = diag([1.0 1.0].^2);

nl = numel(noise);

tail  = t >= tTail;                          % settled part of the response
Tgrid = logspace(log10(2),log10(600),60);    % [s] search range for time constants

parname = {'K11','K12','K21','K22','T1','T2','T3','T4'};
np      = numel(parname);

%% Linearized model at the operating point (reference values)
hs   = xs./(rho*At);
Tlin = (At./ap).*sqrt(2*hs/g);                                        % [s]
Klin = [gamma1*Tlin(1)/(rho*At(1))      (1-gamma2)*Tlin(1)/(rho*At(1));
        (1-gamma1)*Tlin(2)/(rho*At(2))  gamma2*Tlin(2)/(rho*At(2))];  % [cm/(cm^3/s)]
thetaLin = [Klin(1,1); Klin(1,2); Klin(2,1); Klin(2,2); Tlin];

%% Deterministic step responses, identified with the same method
thetaDet = zeros(np,ns);
opts = odeset('RelTol',1e-8,'AbsTol',1e-6);
for j = 1:ns
    Sj = cell(1,nu);
    for i = 1:nu
        ustep    = us;
        du       = stepPct(j)*us(i);
        ustep(i) = us(i) + du;
        [~,X]    = ode15s(@(tt,xx) QuadrupleTankProcess(tt,xx,ustep,ds,p), t, xs, opts);
        Sj{i}    = (X(:,1:2)'./(rho*At(1:2)) - ys)/du;
    end
    thetaDet(:,j) = identifyModel(t,Sj{1},Sj{2},tail,Tgrid);
end

%% Random numbers for all realizations
% W : nw x ((N-1)*nsub+1) x Nmc  Wiener paths
% E : ny x N x Nmc               standardized measurement noise
[W,Tw] = StdWienerProcess(tfinal,(N-1)*nsub,nw,Nmc,seed);
rng(seed+1);
E = randn(ny,N,Nmc);

kSample = 1:nsub:(N-1)*nsub+1;     % indices of the measurement times in Tw

%% Monte Carlo identification
% theta(:,l,j,r): [K11 K12 K21 K22 T1 T2 T3 T4] for noise level l,
%                 step size j and realization r
theta = zeros(np,nl,ns,Nmc);
tic
for r = 1:Nmc
    for l = 1:nl
        psde = [p; noise(l).sigma];              % layout expected by QuadrupleTankDiffusion
        v    = chol(noise(l).Rvv,'lower')*E(:,:,r);
        for j = 1:ns
            Sj = cell(1,nu);
            for i = 1:nu
                ustep    = us;
                du       = stepPct(j)*us(i);
                ustep(i) = us(i) + du;

                % u and d are constant after the step -> one SDE solve
                X = SDEsolverExplicitExplicit(@QuadrupleTankProcess, ...
                                              @QuadrupleTankDiffusion, ...
                                              Tw, xs, W(:,:,r), ustep, ds, psde);
                y     = X(1:2,kSample)./(rho*At(1:2)) + v;   % y = g(x,p) + v
                Sj{i} = (y - ys)/du;                         % normalized response
            end
            theta(:,l,j,r) = identifyModel(t,Sj{1},Sj{2},tail,Tgrid);
        end
    end
end
fprintf('Monte Carlo identification: %d realizations in %.0f s\n',Nmc,toc);

%% Statistics over the realizations
thetaMean = mean(theta,4);                                    % np x nl x ns
thetaStd  = std(theta,0,4);
relRmse   = sqrt(mean((theta - thetaLin).^2,4))./thetaLin;    % relative to linearized model
relDet    = abs(thetaDet - thetaLin)./thetaLin;               % deterministic: bias only

% Fraction of time constant fits that end at the limit of the search range
atLimit = mean(theta(5:8,:,:,:) <= 1.001*Tgrid(1) | theta(5:8,:,:,:) >= 0.999*Tgrid(end),4);

%% Print results
fprintf('\nMean +- standard deviation of the identified parameters (%d realizations)\n',Nmc);
fprintf('Gains in cm/(cm^3/s), time constants in s\n');
fprintf('%-7s %-5s','noise','step');
fprintf(' %17s',parname{:});
fprintf('\n');
fprintf('%-7s %-5s','linear','');
fprintf(' %17.4g',thetaLin);
fprintf('\n');
for j = 1:ns
    fprintf('%-7s +%2.0f%% ','determ.',100*stepPct(j));
    fprintf(' %17.4g',thetaDet(:,j));
    fprintf('\n');
end
for l = 1:nl
    for j = 1:ns
        fprintf('%-7s +%2.0f%% ',noise(l).name,100*stepPct(j));
        for k = 1:np
            fprintf(' %8.4g +- %-6.2g',thetaMean(k,l,j),thetaStd(k,l,j));
        end
        fprintf('\n');
    end
end

fprintf('\nRelative RMS error with respect to the linearized model [%%]\n');
fprintf('%-7s %-5s','noise','step');
fprintf(' %7s',parname{:});
fprintf('\n');
for j = 1:ns
    fprintf('%-7s +%2.0f%% ','determ.',100*stepPct(j));
    fprintf(' %7.1f',100*relDet(:,j));
    fprintf('\n');
end
for l = 1:nl
    for j = 1:ns
        fprintf('%-7s +%2.0f%% ',noise(l).name,100*stepPct(j));
        fprintf(' %7.1f',100*relRmse(:,l,j));
        fprintf('\n');
    end
end

fprintf('\nFraction of fits at the limit of the search range [%g, %g] s [%%]\n',Tgrid(1),Tgrid(end));
fprintf('%-7s %-5s %7s %7s %7s %7s\n','noise','step',parname{5:8});
for l = 1:nl
    for j = 1:ns
        fprintf('%-7s +%2.0f%%  %7.1f %7.1f %7.1f %7.1f\n', ...
                noise(l).name,100*stepPct(j),100*atLimit(:,l,j));
    end
end

%% Figure: relative RMS error against step size, one line per noise level
% Row 1: gains. Row 2: time constants. Logarithmic error axis, common
% limits within each row.
colNoise = [134 182 239;      % low     (light -> dark = low -> high noise)
             42 120 214;      % medium
             13  54 107]/255; % high
colDet   = [0.35 0.35 0.35];
mk       = {'o','s','^'};
ptitle   = {'K_{11}','K_{12}','K_{21}','K_{22}','T_1','T_2','T_3','T_4'};
xs_      = 1:ns;                                  % step sizes as categories
xlab     = arrayfun(@(s) sprintf('+%.0f%%',100*s),stepPct,'UniformOutput',false);

fig = figure('Name','4.5 identification from noisy step responses','Position',[100 100 1100 600]);
tl  = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
hl  = gobjects(1,nl+1);
for k = 1:np
    ax = nexttile(k);
    hl(1) = semilogy(xs_,100*relDet(k,:),'--x','Color',colDet,'LineWidth',1.2,'MarkerSize',7);
    hold on
    for l = 1:nl
        hl(l+1) = semilogy(xs_,100*squeeze(relRmse(k,l,:)),['-' mk{l}], ...
                           'Color',colNoise(l,:),'MarkerFaceColor',colNoise(l,:), ...
                           'LineWidth',1.6,'MarkerSize',6);
    end
    hold off
    grid on; box on
    set(ax,'XTick',xs_,'XTickLabel',xlab,'YMinorGrid','off')
    xlim([0.7 ns+0.3])
    title(ptitle{k})
    if k > 4, xlabel('step size'), end
    if k == 1 || k == 5, ylabel('relative RMS error [%]'), end
end
% Common error axis within each row
for row = 1:2
    idx = (row-1)*4 + (1:4);
    v   = 100*[reshape(relRmse(idx,:,:),1,[]) reshape(relDet(idx,:),1,[])];
    % Lower limit at least 0.1%: smaller errors leave the axis at the bottom.
    % Upper limit: largest error rounded up within its decade.
    dec = 10^floor(log10(max(v)));
    yl  = [max(10^floor(log10(min(v))),0.1) ceil(max(v)/dec)*dec];
    for k = idx
        ylim(nexttile(k),yl)
    end
end
lg = legend(hl,[{'deterministic'} arrayfun(@(n) [n.name ' noise'],noise,'UniformOutput',false)], ...
            'Orientation','horizontal');
lg.Layout.Tile = 'north';

%% Save results
datadir = fullfile(here,'data');
if ~exist(datadir,'dir'), mkdir(datadir), end
save(fullfile(datadir,'noise_identification.mat'), ...
     'stepPct','noise','Nmc','tTail','parname','thetaLin','thetaDet', ...
     'thetaMean','thetaStd','relRmse','relDet','atLimit');

plotdir = fullfile(here,'plots');
if ~exist(plotdir,'dir'), mkdir(plotdir), end
try, theme(fig,'light'), catch, end   % report figures on white (R2025a+)
exportgraphics(fig,fullfile(plotdir,'sim_4_5_noise_identification.png'),'Resolution',150)
exportgraphics(fig,fullfile(plotdir,'sim_4_5_noise_identification.pdf'),'ContentType','vector')

%% Local functions
function theta = identifyModel(t,SF1,SF2,tail,Tgrid)
% Identify [K11 K12 K21 K22 T1 T2 T3 T4] from the normalized responses
% SF1 (2 x N: h1, h2 for a step in F1) and SF2 (2 x N for a step in F2).
S11 = SF1(1,:);   S21 = SF1(2,:);
S12 = SF2(1,:);   S22 = SF2(2,:);

% 1) gains: average over the settled part of the response
K11 = mean(S11(tail));   K21 = mean(S21(tail));
K12 = mean(S12(tail));   K22 = mean(S22(tail));

% 2) T1, T2: first-order fit to the direct responses
T1 = fitTimeConstant(@(T) sum((S11 - firstOrderStep(t,K11,T)).^2), Tgrid);
T2 = fitTimeConstant(@(T) sum((S22 - firstOrderStep(t,K22,T)).^2), Tgrid);

% 3) T3, T4: second-order fit to the cross responses
T3 = fitTimeConstant(@(T) sum((S12 - secondOrderStep(t,K12,T1,T)).^2), Tgrid);
T4 = fitTimeConstant(@(T) sum((S21 - secondOrderStep(t,K21,T2,T)).^2), Tgrid);

theta = [K11; K12; K21; K22; T1; T2; T3; T4];
end

function T = fitTimeConstant(sse,Tgrid)
% Least-squares estimate of one time constant: coarse search over Tgrid,
% then refinement with fminbnd between the neighbours of the best point.
% With noise the sum of squares can have several local minima, so a
% purely local search is not reliable.
v      = arrayfun(sse,Tgrid);
[~,k]  = min(v);
lo     = Tgrid(max(k-1,1));
hi     = Tgrid(min(k+1,numel(Tgrid)));
T      = fminbnd(sse,lo,hi,optimset('TolX',1e-3));
end

function y = firstOrderStep(t,K,T)
% Step response of K/(T s+1)
y = K*(1 - exp(-t/T));
end

function y = secondOrderStep(t,K,Ta,Tb)
% Step response of K/((Ta s+1)(Tb s+1))
if abs(Ta - Tb) < 1e-8*max(Ta,1)
    y = K*(1 - (1 + t/Ta).*exp(-t/Ta));
else
    y = K*(1 - (Ta*exp(-t/Ta) - Tb*exp(-t/Tb))/(Ta - Tb));
end
end
