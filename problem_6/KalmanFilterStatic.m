function [xf,xp,e] = KalmanFilterStatic(kf,xp,y,u)
% One step of the static (stationary) Kalman filter: same recursion as
% KalmanFilterDynamic with the constant gain kf.Kfx from KalmanFilterDesign.
% In:  xp    one-step prediction x_{k|k-1},  y, u  measurement and input
% Out: xf    filtered x_{k|k}
%      xp    prediction x_{k+1|k} for the next call
%      e     innovation

% Filtering
e  = y - kf.C*xp;
xf = xp + kf.Kfx*e;

% One-step prediction
if nargout > 1
    xp = kf.A*xf + kf.B*u;
end
end
