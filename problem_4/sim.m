
stepPct = [0.1 0.25 0.5];
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
