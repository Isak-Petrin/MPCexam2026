function [u,ctrl] = PIDController(ctrl,r,y,Ts)
% Discrete-time decentralized PID: one SISO loop per element of y and u.
%
%   u_i = ubar_i + Kc_i * ( e_i + (1/tauI_i) * sum(e_i) Ts - tauD_i * dy_i/dt )
%
%   P   : tauI = Inf, tauD = 0
%   PI  : tauD = 0
%   PID : Kc, tauI, tauD all active
%


    if ~isfield(ctrl,'I') || isempty(ctrl.I)
        ctrl.I     = zeros(size(y));
        ctrl.D     = zeros(size(y));
        ctrl.yprev = y;
    end

    e = r - y;

    % Proportional
    P = ctrl.Kc.*e;

    % Filtered derivative on the measurement (zero when tauD = 0)
    tauf = ctrl.tauD./ctrl.N;
    a    = tauf./(tauf + Ts);
    D    = a.*ctrl.D - ctrl.Kc.*ctrl.tauD./(tauf + Ts).*(y - ctrl.yprev);

    % Integral (zero increment when tauI = Inf)
    I = ctrl.I + ctrl.Kc.*(Ts./ctrl.tauI).*e;

    % Output with saturation
    v = ctrl.ubar + P + I + D;
    u = min(max(v,ctrl.umin),ctrl.umax);

    % Anti-windup: keep the previous integral in saturated loops
    sat    = (u ~= v);
    I(sat) = ctrl.I(sat);
    v = ctrl.ubar + P + I + D;
    u = min(max(v,ctrl.umin),ctrl.umax);

    ctrl.I     = I;
    ctrl.D     = D;
    ctrl.yprev = y;
end
