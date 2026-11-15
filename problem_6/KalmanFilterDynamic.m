function [xf,Pf,xp,Pp,e,Re] = KalmanFilterDynamic(kf,xp,Pp,y,u)
% One step of the dynamic (time-varying) Kalman filter.
% In:  xp, Pp   one-step prediction x_{k|k-1}, P_{k|k-1}
%      y, u     measurement y_k and input u_k
% Out: xf, Pf   filtered x_{k|k}, P_{k|k}
%      xp, Pp   prediction x_{k+1|k}, P_{k+1|k} for the next call
%      e, Re    innovation and its covariance
% In closed loop u_k depends on xf: call [xf,Pf] = KalmanFilterDynamic(kf,xp,Pp,y),
% compute u_k, then [xp,Pp] = KalmanFilterPredict(kf,xf,Pf,u).

C = kf.C;

% Filtering (measurement update)
e   = y - C*xp;
Re  = C*Pp*C' + kf.R;
Kfx = Pp*C'/Re;
xf  = xp + Kfx*e;
IKC = eye(length(xp)) - Kfx*C;
Pf  = IKC*Pp*IKC' + Kfx*kf.R*Kfx';      % Joseph form, keeps Pf symmetric

% One-step prediction (time update)
if nargout > 2
    xp = kf.A*xf + kf.B*u;
    Pp = kf.A*Pf*kf.A' + kf.Qw;
end
end
