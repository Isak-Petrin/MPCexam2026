%% Problem 4.2 - stochastic step responses (process + measurement noise)
% Steps of +10%, +25% and +50% in F1 and in F2 (one input at a time) on
% the stochastic nonlinear model (SDE, model 2.3) with measurement noise:
%
%   dx = f(x,u,d,p) dt + sigma(x,u,d,p) dw      process noise on F3, F4
%   y  = g(x,p) + v,      v ~ N(0,Rvv)           measurement noise
%
% The drift uses the nominal disturbance d = ds, so the simulations start
% in the same steady state xs as Problem 4.1 and the normalized step
% responses are directly comparable with the deterministic ones.
%
% Three noise levels (low / medium / high). The medium level equals the
% noise parameters used in Problem 2 (sigma = 20, std(v) = 0.5 cm).
%
% Common random numbers: within one noise level the SAME Wiener path and
% the SAME measurement noise sequence are used for all step sizes and
% both inputs, so differences between curves are caused by the step size
% only. Across noise levels the same standardized noise is scaled.
%
% Normalized step responses (Problem 4.3)
%     S(t) = (y(t) - y_s) / (u - u_s)        [cm / (cm^3/s)]
% are plotted and everything is saved to data/steps_stochastic.mat.
%
% Layout follows Azam (Fig. 2.6-2.8): one figure pair per noise level.

clc
clear
close all

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','problem_2','problem_2_1'));
addpath(fullfile(here,'..','problem_2','problem_2_4'));

QuadrupleTankParameters   % p, us, ds, xs, ys, rho, At, seed, ...

%% Settings
stepPct = [0.10 0.25 0.50];
Ts      = 10;            % [s] sampling time (measurements)
tfinal  = 20*60;         % [s]
t       = 0:Ts:tfinal;
tmin    = t/60;
N       = length(t);
nsub    = 10;            % Euler-Maruyama sub-steps per sample (dt = 1 s)

nx = 4;
nu = 2;
ny = 2;
nw = 2;
ns = length(stepPct);

%% Noise levels  (CHOICE)
noise(1).name  = 'low';
noise(1).sigma = [ 5;  5];            % [cm^3/s^0.5] diffusion on F3, F4
noise(1).Rvv   = diag([0.1 0.1].^2);  % [cm^2]

noise(2).name  = 'medium';            % = Problem 2 values
noise(2).sigma = [20; 20];
noise(2).Rvv   = diag([0.5 0.5].^2);

noise(3).name  = 'high';
noise(3).sigma = [50; 50];
noise(3).Rvv   = diag([1.0 1.0].^2);

nl = numel(noise);

%% Common random numbers (one realization, reused everywhere)
[W,Tw,~] = StdWienerProcess(tfinal,(N-1)*nsub,nw,1,seed);   % nw x ((N-1)*nsub+1)
rng(seed+1);
e = randn(ny,N);                                            % standardized meas. noise

%% Deterministic reference (same as Problem 4.1, measured outputs only)
Sdet = cell(nu,ns);
du   = zeros(nu,ns);
opts = odeset('RelTol',1e-8,'AbsTol',1e-6);
for i = 1:nu
    for j = 1:ns
        ustep    = us;
        du(i,j)  = stepPct(j)*us(i);
        ustep(i) = us(i) + du(i,j);
        [~,X] = ode15s(@(tt,xx) QuadrupleTankProcess(tt,xx,ustep,ds,p), t, xs, opts);
        Sdet{i,j} = (X(:,1:2)'./(rho*At(1:2)) - ys)/du(i,j);
    end
end

%% Stochastic simulations
% Y{l,i,j} : 2 x N measurements y for noise level l, input i, step j
% Z{l,i,j} : 2 x N true levels   z (noise-free outputs of the noisy process)
% S{l,i,j} : 2 x N normalized measured step response
Y = cell(nl,nu,ns);
Z = cell(nl,nu,ns);
S = cell(nl,nu,ns);

for l = 1:nl
    psde = [p; noise(l).sigma];       % layout expected by QuadrupleTankDiffusion
    v    = chol(noise(l).Rvv,'lower')*e;

    for i = 1:nu
        for j = 1:ns
            ustep    = us;
            ustep(i) = us(i) + du(i,j);

            x = zeros(nx,N);
            x(:,1) = xs;
            for k = 1:N-1
                idx = (k-1)*nsub + (1:nsub+1);
                Xk  = SDEsolverExplicitExplicit(@QuadrupleTankProcess, ...
                                                @QuadrupleTankDiffusion, ...
                                                Tw(idx), x(:,k), W(:,idx), ...
                                                ustep, ds, psde);
                x(:,k+1) = Xk(:,end);
            end

            z = zeros(ny,N);
            for k = 1:N
                z(:,k) = QuadrupleTankOutput(x(:,k),p);
            end
            y = z + v;                                 % y = g(x,p) + v

            Z{l,i,j} = z;
            Y{l,i,j} = y;
            S{l,i,j} = (y - ys)/du(i,j);
        end
    end
end

%% Print: noise in the normalized responses (std of S - Sdet over the last 10 min)
fprintf('\nStd of (S_noisy - S_det) over t >= 10 min  [cm/(cm^3/s)]\n');
fprintf('%-8s %-6s %-6s %10s %10s\n','noise','input','step','y1','y2');
late = t >= 10*60;
for l = 1:nl
    for i = 1:nu
        for j = 1:ns
            r = S{l,i,j}(:,late) - Sdet{i,j}(:,late);
            fprintf('%-8s F%-5d %3.0f%%   %10.4f %10.4f\n', ...
                    noise(l).name,i,100*stepPct(j),std(r,0,2));
        end
    end
end

%% Plot settings
cols    = [0.00 0.45 0.74;
           0.85 0.33 0.10;
           0.47 0.67 0.19];
legtxt  = arrayfun(@(s) sprintf('+%.0f%%',100*s),stepPct,'UniformOutput',false);
outname = {'y_1 = h_1','y_2 = h_2'};
outshort = {'y_1','y_2'};

plotdir = fullfile(here,'plots');
if ~exist(plotdir,'dir'), mkdir(plotdir), end

for l = 1:nl
    %% Step responses (measured levels)
    f1 = figure('Name',sprintf('4.2 step responses, %s noise',noise(l).name), ...
                'Position',[100 100 900 550]);
    for i = 1:nu
        for r = 1:ny
            subplot(ny,nu,(r-1)*nu+i)
            hold on
            for j = 1:ns
                plot(tmin,Y{l,i,j}(r,:),'Color',cols(j,:),'LineWidth',1.0)
            end
            plot(tmin([1 end]),ys(r)*[1 1],'k--','LineWidth',0.8)
            hold off
            grid on; box on
            xlim(tmin([1 end]))
            ylabel([outname{r} ' [cm]'])
            if r == 1
                title(sprintf('Step in F_%d (%s noise)',i,noise(l).name))
            end
            if r == ny, xlabel('time [min]'), end
            if r == 1 && i == nu
                legend([legtxt {'steady state'}],'Location','southeast')
            end
        end
    end

    %% Normalized step responses (with deterministic 10% as reference)
    f2 = figure('Name',sprintf('4.2 normalized step responses, %s noise',noise(l).name), ...
                'Position',[150 100 900 550]);
    for i = 1:nu
        for r = 1:ny
            subplot(ny,nu,(r-1)*nu+i)
            hold on
            for j = 1:ns
                plot(tmin,S{l,i,j}(r,:),'Color',cols(j,:),'LineWidth',1.0)
            end
            plot(tmin,Sdet{i,1}(r,:),'k--','LineWidth',1.2)
            hold off
            grid on; box on
            xlim(tmin([1 end]))
            ylabel(sprintf('\\Delta%s/\\DeltaF_%d',outshort{r},i))
            if r == 1
                title(sprintf('Normalized step in F_%d (%s noise)',i,noise(l).name))
            end
            if r == ny, xlabel('time [min]'), end
            if r == 1 && i == nu
                legend([legtxt {'deterministic'}],'Location','southeast')
            end
        end
    end

    %% Save figures
    figs  = [f1 f2];
    names = {sprintf('sim_4_2_steps_%s_noise',noise(l).name), ...
             sprintf('sim_4_2_steps_%s_noise_normalized',noise(l).name)};
    for k = 1:2
        try, theme(figs(k),'light'), catch, end
        exportgraphics(figs(k),fullfile(plotdir,[names{k} '.png']),'Resolution',150)
        exportgraphics(figs(k),fullfile(plotdir,[names{k} '.pdf']),'ContentType','vector')
    end
end

%% Save data for Problem 4.3 / 4.4
datadir = fullfile(here,'data');
if ~exist(datadir,'dir'), mkdir(datadir), end
save(fullfile(datadir,'steps_stochastic.mat'), ...
     't','tmin','Ts','stepPct','us','ds','xs','ys','du','noise','Y','Z','S','Sdet');
