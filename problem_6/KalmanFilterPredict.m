function [xp,Pp] = KalmanFilterPredict(kf,xf,Pf,u)
% One-step prediction (time update) of the Kalman filter:
%   x_{k+1|k} = A x_{k|k} + B u_k,   P_{k+1|k} = A P_{k|k} A' + G Q G'
% Use it after KalmanFilterDynamic(kf,xp,Pp,y) when u_k is computed from
% the filtered estimate (closed loop with a controller).

xp = kf.A*xf + kf.B*u;
if nargout > 1
    Pp = kf.A*Pf*kf.A' + kf.Qw;
end
end
