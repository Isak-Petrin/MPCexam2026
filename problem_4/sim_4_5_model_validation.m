%% Problem 4.5 - validation of the identified transfer function model
% Compares the identified model
%
%            [ K11/(T1 s+1)               K12/((T1 s+1)(T3 s+1)) ]
%     G(s) = [                                                    ]
%            [ K21/((T2 s+1)(T4 s+1))     K22/(T2 s+1)            ]
%
% with the deterministic normalized step responses it was identified from,
% in all four channels, and reports the fit error per channel.
%
% The identification is the same as in sim_4_4_curvefit_T3_T4.m:
%   1) gains  K_ij = S_ij(t_end)                         (end value)
%   2) T1, T2      = time where S11, S22 reach 63.2% of K11, K22
%   3) T3, T4      = least-squares fit of the second-order step responses
% Step 3 solves the same least-squares problems as sim_4_4 (lsqcurvefit),
% but with fminbnd, and the model step responses are evaluated from their
% analytical expressions. The script therefore runs in base MATLAB
% (no Optimization Toolbox, no Control System Toolbox).
%
% Output:
%   - printed table of K, T, RMSE and maximum error for the +10/+25/+50%
%     steps (numbers for the accuracy discussion)
%   - plots/sim_4_5_model_validation.{png,pdf}: model against simulation
%     and residuals for the +10% step (the identified model)

clc
clear
close all

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','problem_2','problem_2_1'));

QuadrupleTankParameters   % p, us, ds, xs, rho, At, ...

%% Settings (same grid as Problem 4.1)
stepPct = [0.10 0.25 0.50];
Ts      = 10;                 % [s]
tfinal  = 20*60;              % [s]
t       = 0:Ts:tfinal;
tmin    = t/60;
nu      = 2;
ny      = 2;
ns      = length(stepPct);
jId     = 1;                  % step used for the identified model (+10%)

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

%% Identification and model step responses for every step size
% Index (r,i,j) = (output h_r, input F_i, step size j)
K    = zeros(ny,nu,ns);       % gains
T    = zeros(4,ns);           % [T1; T2; T3; T4]
Sd   = cell(ny,nu,ns);        % simulated normalized step response (data)
Sm   = cell(ny,nu,ns);        % model step response
rmse = zeros(ny,nu,ns);       % root mean square error
emax = zeros(ny,nu,ns);       % maximum |error| relative to the gain

fbopts = optimset('TolX',1e-8);
for j = 1:ns
    for i = 1:nu
        for r = 1:ny
            Sd{r,i,j} = S{i,j}(r,:);
            K(r,i,j)  = Sd{r,i,j}(end);              % 1) gain = end value
        end
    end

    % 2) T1, T2 from the 63.2% rule
    T1 = time63(t,Sd{1,1,j});
    T2 = time63(t,Sd{2,2,j});

    % 3) T3, T4 by least squares (K and T1, T2 fixed)
    T3 = fminbnd(@(x) sum((Sd{1,2,j} - secondOrderStep(t,K(1,2,j),T1,x)).^2), 1, 300, fbopts);
    T4 = fminbnd(@(x) sum((Sd{2,1,j} - secondOrderStep(t,K(2,1,j),T2,x)).^2), 1, 300, fbopts);
    T(:,j) = [T1; T2; T3; T4];

    % Model step responses of the four channels
    Sm{1,1,j} = firstOrderStep(t,K(1,1,j),T1);
    Sm{1,2,j} = secondOrderStep(t,K(1,2,j),T1,T3);
    Sm{2,1,j} = secondOrderStep(t,K(2,1,j),T2,T4);
    Sm{2,2,j} = firstOrderStep(t,K(2,2,j),T2);

    % Fit error per channel
    for i = 1:nu
        for r = 1:ny
            e           = Sd{r,i,j} - Sm{r,i,j};
            rmse(r,i,j) = sqrt(mean(e.^2));
            emax(r,i,j) = max(abs(e))/K(r,i,j);
        end
    end
end

%% Print results
fprintf('\nIdentified model per step size (deterministic step responses)\n');
fprintf('%-6s %8s %8s %8s %8s | %7s %7s %7s %7s\n', ...
        'step','K11','K12','K21','K22','T1','T2','T3','T4');
for j = 1:ns
    fprintf('+%3.0f%%  %8.4f %8.4f %8.4f %8.4f | %7.1f %7.1f %7.1f %7.1f\n', ...
            100*stepPct(j),K(1,1,j),K(1,2,j),K(2,1,j),K(2,2,j),T(:,j));
end

fprintf('\nFit error of the model identified from each step size\n');
fprintf('%-6s | %9s %9s %9s %9s | %7s %7s %7s %7s\n', ...
        'step','RMSE S11','RMSE S12','RMSE S21','RMSE S22', ...
        'max S11','max S12','max S21','max S22');
for j = 1:ns
    fprintf('+%3.0f%%  | %9.2e %9.2e %9.2e %9.2e | %6.2f%% %6.2f%% %6.2f%% %6.2f%%\n', ...
            100*stepPct(j),rmse(1,1,j),rmse(1,2,j),rmse(2,1,j),rmse(2,2,j), ...
            100*emax(1,1,j),100*emax(1,2,j),100*emax(2,1,j),100*emax(2,2,j));
end
fprintf('(max = largest |simulation - model| in percent of the channel gain)\n');

%% Figure: identified model against simulation, all four channels (+10% step)
% Columns: input F_1, F_2. Rows: output h_1, h_2. Under each response the
% residual (simulation - model) is shown on a common scale.
cols     = [0.00 0.45 0.74; 0.85 0.33 0.10];
resScale = 1e5;
resMax   = 0;
for i = 1:nu
    for r = 1:ny
        resMax = max(resMax, max(abs(Sd{r,i,jId} - Sm{r,i,jId})));
    end
end
resLim = 5*ceil(resMax*resScale/5);       % common residual axis limit

% Parameters shown in each panel
Tid  = T(:,jId);
ptxt = cell(ny,nu);
ptxt{1,1} = sprintf('K_{11} = %.4f\nT_1 = %.1f s',               K(1,1,jId),Tid(1));
ptxt{1,2} = sprintf('K_{12} = %.4f\nT_1 = %.1f s, T_3 = %.1f s', K(1,2,jId),Tid(1),Tid(3));
ptxt{2,1} = sprintf('K_{21} = %.4f\nT_2 = %.1f s, T_4 = %.1f s', K(2,1,jId),Tid(2),Tid(4));
ptxt{2,2} = sprintf('K_{22} = %.4f\nT_2 = %.1f s',               K(2,2,jId),Tid(2));

fig = figure('Name','4.5 identified model vs simulation','Position',[100 100 950 780]);
tl  = tiledlayout(6,2,'TileSpacing','compact','Padding','compact');
for r = 1:ny
    for i = 1:nu
        y  = Sd{r,i,jId};
        ym = Sm{r,i,jId};

        % Response: rows 1-2 (h_1) and 4-5 (h_2) of the 6 x 2 layout
        nexttile((r-1)*6 + i,[2 1])
        hd = plot(tmin,y,'o','Color',cols(1,:),'MarkerSize',4); hold on
        hm = plot(tmin,ym,'-','Color',cols(2,:),'LineWidth',1.4); hold off
        grid on; box on
        xlim(tmin([1 end])); ylim([0 1.08*K(r,i,jId)])
        set(gca,'XTickLabel',[])
        ylabel(sprintf('\\Deltah_%d/\\DeltaF_%d',r,i))
        if r == 1, title(sprintf('Step in F_%d',i)), end
        text(0.97,0.08,sprintf('%s\nRMSE = %.1f\\cdot10^{-5}',ptxt{r,i},rmse(r,i,jId)*resScale), ...
             'Units','normalized','HorizontalAlignment','right', ...
             'VerticalAlignment','bottom','FontSize',9)
        if r == 1 && i == 1
            lg = legend([hd hm],{'nonlinear simulation (+10%)','identified model'}, ...
                        'Orientation','horizontal');
            lg.Layout.Tile = 'north';
        end

        % Residual: row 3 (h_1) and row 6 (h_2)
        nexttile((r-1)*6 + 4 + i)
        plot(tmin,(y - ym)*resScale,'k-','LineWidth',1.0)
        grid on; box on
        xlim(tmin([1 end])); ylim(resLim*[-1 1])
        ylabel({'residual','[10^{-5}]'})
        if r == ny
            xlabel('time [min]')
        else
            set(gca,'XTickLabel',[])
        end
    end
end

%% Save figure
plotdir = fullfile(here,'plots');
if ~exist(plotdir,'dir'), mkdir(plotdir), end
try, theme(fig,'light'), catch, end   % report figures on white (R2025a+)
exportgraphics(fig,fullfile(plotdir,'sim_4_5_model_validation.png'),'Resolution',150)
exportgraphics(fig,fullfile(plotdir,'sim_4_5_model_validation.pdf'),'ContentType','vector')

%% Local functions
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

function T = time63(t,s)
% Time at which s first reaches 63.2% of its final value
target = 0.632*s(end);
k = find(s >= target, 1, 'first');
T = t(k-1) + (target - s(k-1))*(t(k) - t(k-1))/(s(k) - s(k-1));
end
