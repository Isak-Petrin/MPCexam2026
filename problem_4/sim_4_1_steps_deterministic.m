%% Problem 4.1 - deterministic step responses in the manipulated variables
% Steps of +10%, +25% and +50% in F1 and in F2

clc
clear
close all

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','problem_2','problem_2_1'));

QuadrupleTankParameters   % p, us, ds, xs, zs, ...

%% Settings
stepPct = [0.10 0.25 0.50];
Ts      = 10;         %  same grid as Prob 2
tfinal  = 20*60;      
t       = 0:Ts:tfinal;
tmin    = t/60;
N       = length(t);

nx = 4;
nu = 2;
nz = 2;

%% Step responses: xNL{i,j}, zNL{i,j} = step of stepPct(j) in input i
xNL = cell(nu,length(stepPct));
zNL = cell(nu,length(stepPct));

for i = 1:nu
    for j = 1:length(stepPct)
        ustep    = us;
        ustep(i) = us(i)*(1+stepPct(j));

        x = zeros(nx,N);
        z = zeros(nz,N);
        x(:,1) = xs;

        for k = 1:N-1
            z(:,k) = QuadrupleTankOutput(x(:,k),p);

            [~,Xk] = ode15s(@QuadrupleTankProcess, ...
                            [t(k) t(k+1)], ...
                            x(:,k), ...
                            [], ...
                            ustep, ds, p);

            x(:,k+1) = Xk(end,:)';
        end

        z(:,N) = QuadrupleTankOutput(x(:,N),p);

        xNL{i,j} = x;
        zNL{i,j} = z;
    end
end


for i = 1:nu
    figure('Name',sprintf('4.1 steps in F_%d',i));
    for j = 1:length(stepPct)
        subplot(3,1,j)
        plot(tmin,zNL{i,j}(1,:),'LineWidth',1.2)
        hold on
        plot(tmin,zNL{i,j}(2,:),'LineWidth',1.2)
        yline(zs(1),'--')
        yline(zs(2),'--')
        xlabel('time [min]')
        ylabel('level [cm]')
        title(sprintf('Nonlinear response: +%.0f%% step in F_%d',100*stepPct(j),i))
        legend('h_1','h_2','h_{1,s}','h_{2,s}','Location','best')
        grid on
    end
end

plotdir = fullfile(here,'plots');
if ~exist(plotdir,'dir')
    mkdir(plotdir)
end
for i = 1:nu
    f = figure(i);
    try, theme(f,'light'), catch, end   % report figures on white (R2025a+)
    exportgraphics(f,fullfile(plotdir,sprintf('sim_4_1_steps_F%d.png',i)),'Resolution',150)
    exportgraphics(f,fullfile(plotdir,sprintf('sim_4_1_steps_F%d.pdf',i)),'ContentType','vector')
end
