%% Problem 4.1 - deterministic step responses in the manipulated variables
% Steps of +10%, +25% and +50% in F1 and in F2 (one input at a time) on
% the deterministic nonlinear model (2.1), starting from the steady state
% (us, ds, xs) defined in QuadrupleTankParameters.
%
% Also computes the normalized step responses (Problem 4.3)
%     S(t) = (h(t) - h_s) / (u - u_s)        [cm / (cm^3/s)]
% and saves everything to data/steps_deterministic.mat for Problem 4.4.
%
% Layout of the figures follows Azam (Fig. 2.5): one column per input,
% one row per tank level, the three step sizes overlaid in each subplot.

clc
clear
close all

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','problem_2','problem_2_1'));

QuadrupleTankParameters   % p, us, ds, xs, zs, rho, At, ...

%% Settings
stepPct = [0.10 0.25 0.50];   % relative step sizes
Ts      = 10;                 % [s] output grid, same as Problem 2
tfinal  = 20*60;              % [s] long enough to reach the new steady state
t       = 0:Ts:tfinal;
tmin    = t/60;
N       = length(t);

nx = 4;
nu = 2;
ns = length(stepPct);

hs = xs./(rho*At);            % steady-state levels of all four tanks [cm]

%% Simulation
% H{i,j}  : 4 x N tank levels for a step of stepPct(j) in input i
% S{i,j}  : 4 x N normalized step response
% du(i,j) : absolute step size [cm^3/s]
H  = cell(nu,ns);
S  = cell(nu,ns);
du = zeros(nu,ns);

opts = odeset('RelTol',1e-8,'AbsTol',1e-6);

for i = 1:nu
    for j = 1:ns
        ustep    = us;
        du(i,j)  = stepPct(j)*us(i);
        ustep(i) = us(i) + du(i,j);

        % u and d are constant after the step -> one ODE solve on the grid t
        [~,X] = ode15s(@(tt,xx) QuadrupleTankProcess(tt,xx,ustep,ds,p), t, xs, opts);
        X = X';                                   % 4 x N

        H{i,j} = X./(rho*At);                     % levels [cm]
        S{i,j} = (H{i,j} - hs)/du(i,j);           % normalized step response
        S{i,j}(abs(S{i,j}) < 1e-10) = 0;          % remove round-off (h3 for F1, h4 for F2)
    end
end

%% Print steady-state gains (end value of the normalized responses)
fprintf('\nSteady-state gains  K = dh/dF  [cm/(cm^3/s)] at t = %g min\n', tfinal/60);
fprintf('%-8s %-6s %10s %10s %10s %10s\n','input','step','h1','h2','h3','h4');
for i = 1:nu
    for j = 1:ns
        fprintf('F%-7d %3.0f%%   %10.4f %10.4f %10.4f %10.4f\n', ...
                i,100*stepPct(j),S{i,j}(:,end));
    end
end

%% Plot settings
cols    = [0.00 0.45 0.74;    % 10 %
           0.85 0.33 0.10;    % 25 %
           0.47 0.67 0.19];   % 50 %
legtxt  = arrayfun(@(s) sprintf('+%.0f%%',100*s),stepPct,'UniformOutput',false);
outname = {'h_1','h_2','h_3','h_4'};

%% Figure 1: step responses (levels)
fig1 = figure('Name','4.1 deterministic step responses','Position',[100 100 900 900]);
for i = 1:nu
    for r = 1:nx
        subplot(nx,nu,(r-1)*nu+i)
        hold on
        for j = 1:ns
            plot(tmin,H{i,j}(r,:),'Color',cols(j,:),'LineWidth',1.3)
        end
        plot(tmin([1 end]),hs(r)*[1 1],'k--','LineWidth',0.8)
        hold off
        grid on; box on
        xlim(tmin([1 end]))
        ylabel([outname{r} ' [cm]'])
        if r == 1, title(sprintf('Step in F_%d',i)), end
        if r == nx, xlabel('time [min]'), end
        if r == 1 && i == nu
            legend([legtxt {'steady state'}],'Location','southeast')
        end
    end
end

%% Figure 2: normalized step responses
fig2 = figure('Name','4.1 deterministic normalized step responses','Position',[150 100 900 900]);
for i = 1:nu
    for r = 1:nx
        subplot(nx,nu,(r-1)*nu+i)
        hold on
        for j = 1:ns
            plot(tmin,S{i,j}(r,:),'Color',cols(j,:),'LineWidth',1.3)
        end
        hold off
        grid on; box on
        xlim(tmin([1 end]))
        ylabel(sprintf('\\Delta%s/\\DeltaF_%d', outname{r}, i))
        if r == 1, title(sprintf('Normalized step in F_%d',i)), end
        if r == nx, xlabel('time [min]'), end
        if r == 1 && i == nu
            legend(legtxt,'Location','southeast')
        end
    end
end

%% Save data for Problem 4.3 / 4.4
datadir = fullfile(here,'data');
if ~exist(datadir,'dir'), mkdir(datadir), end
save(fullfile(datadir,'steps_deterministic.mat'), ...
     't','tmin','Ts','stepPct','us','ds','xs','hs','du','H','S');

%% Save figures
plotdir = fullfile(here,'plots');
if ~exist(plotdir,'dir'), mkdir(plotdir), end

figs  = [fig1 fig2];
names = {'sim_4_1_steps_deterministic','sim_4_1_steps_deterministic_normalized'};
for k = 1:numel(figs)
    try, theme(figs(k),'light'), catch, end   % report figures on white (R2025a+)
    exportgraphics(figs(k),fullfile(plotdir,[names{k} '.png']),'Resolution',150)
    exportgraphics(figs(k),fullfile(plotdir,[names{k} '.pdf']),'ContentType','vector')
end
