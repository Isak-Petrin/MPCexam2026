%% Problem 5.5 - Markov parameters of the linearized discrete-time model
% Uses the linearized four-tank model from Section 5.1 (Eq. 5.16) and the
% exact zero-order-hold discretization from Section 5.4 (Eq. 5.32).
% Markov parameters are computed the same way as in Problem 4.6:
%   H0 = D,   Hk = C*Ad^(k-1)*Bd,   k = 1,...,N_H-1
% The disturbance matrix E is not part of the Markov parameters for the input
% u; it only enters the disturbance channel (see the commented block below).

clc;
format long g

%% Settings

Ts  = 10;      % Sample time [s]
N_H = 101;     % Number of Markov parameter matrices (H0 to H100)

%% Linearized continuous-time model (Section 5.1)

% Time constants from Table 5.1 [s]
T1 = 95.2;
T2 = 105.5;
T3 = 48.9;
T4 = 51.5;

% Continuous-time matrices (Eq. 5.14 and 5.16), exact form in the time constants
Ac = [-1/T1  0     1/T3  0;
       0    -1/T2  0     1/T4;
       0     0    -1/T3  0;
       0     0     0    -1/T4];

Bc = [0.60  0;
      0     0.70;
      0     0.30;
      0.40  0];

% Sensor matrix (Eq. 5.15 / 5.16): all four tank levels are measured
Cc = 0.002631 * eye(4);

% No direct feedthrough from u to y
Dc = zeros(4,2);

% Disturbance matrix (Eq. 5.14), only needed for the disturbance channel
Ec = [0 0;
      0 0;
      1 0;
      0 1];

%% Exact zero-order-hold discretization (Eq. 5.32, block matrix exponential)

n = size(Ac,1);
m = size(Bc,2);

M  = [Ac, Bc;
      zeros(m,n), zeros(m,m)];
Md = expm(M*Ts);

Ad = Md(1:n,1:n);
Bd = Md(1:n,n+1:n+m);
Cd = Cc;
Dd = Dc;

fprintf('\n============================================\n');
fprintf('PROBLEM 5.5 - DISCRETE-TIME MODEL (Ts = %.1f s)\n', Ts);
fprintf('============================================\n');
disp('A_d ='); disp(Ad);
disp('B_d ='); disp(Bd);
disp('C_d ='); disp(Cd);
disp('D_d ='); disp(Dd);

%% Markov parameters H_k (all four tank levels)

% H is 4 x 2 x N_H: rows = tank levels 1..4, columns = pumps 1 and 2
H55 = zeros(4,2,N_H);

% H0 = D
H55(:,:,1) = Dd;

% Hk = C*Ad^(k-1)*Bd, k = 1,...,N_H-1
v = Bd;
for k = 1:N_H-1
    H55(:,:,k+1) = Cd*v;
    v = Ad*v;
end

%% Lower-tank Markov parameters (tanks 1 and 2), comparable to Problem 4.6

% The same 2x2 structure as the identified model in Problem 4.6
H55_low = H55(1:2,:,:);

fprintf('\n============================================\n');
fprintf('MARKOV PARAMETERS - LOWER TANKS (1 and 2)\n');
fprintf('============================================\n');

for k = 0:5
    fprintf('\nH_%d =\n', k);
    disp(H55_low(:,:,k+1));
end

fprintf('\nH_30 =\n'); disp(H55_low(:,:,31));
fprintf('\nH_100 =\n'); disp(H55_low(:,:,101));

%% Comparison with Problem 4.6 (if its H is in the MATLAB session)

if exist('H','var') && isequal(size(H),[2 2 N_H])
    fprintf('\n============================================\n');
    fprintf('COMPARISON WITH PROBLEM 4.6 (relative difference)\n');
    fprintf('============================================\n');

    for k = [1 2 3 4 5 30 100]
        Href = H(:,:,k+1);          % Problem 4.6 stores H_k at index k+1
        Hnew = H55_low(:,:,k+1);
        relErr = abs(Hnew - Href) ./ abs(Href);
        fprintf('k = %3d: max relative difference = %.2f %%\n', ...
            k, 100*max(relErr(:)));
    end
else
    fprintf('\nProblem 4.6 variable H (2x2x%d) not found in this session.\n', N_H);
    fprintf('Run markov_4_6.m first in the same session to compare.\n');
end

%% Optional: disturbance Markov parameters F_k = C*Ad^(k-1)*Ed (not part of 5.5)
%
% Ed = expm-block of [Ac Ec; 0 0]*Ts (last block column), as in Eq. 5.32:
% Mf = expm([Ac, Ec; zeros(2,n), zeros(2,2)]*Ts);
% Ed = Mf(1:n, n+1:n+2);
% F = zeros(4,2,N_H); v = Ed;
% for k = 1:N_H-1, F(:,:,k+1) = Cd*v; v = Ad*v; end

fprintf('\n============================================\n');
fprintf('PROBLEM 5.5 COMPLETED\n');
fprintf('============================================\n');
