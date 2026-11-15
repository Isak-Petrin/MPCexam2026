
%% Problem 4.6 - Observer canonical realization and discretization
% Uses the identified +10% model from Problem 4.4.
% Discretization using exact zero-order hold (block matrix exponential).

clc;
format long g

%% Settings

Ts = 10;  % Sample time [s]

%% Retrieve identified model from Problem 4.4

if ~exist('Kid','var') || ~exist('Tid','var')
    error(['Kid and Tid are missing. Run the Problem 4.4 ', ...
           'identification script first, in the same MATLAB session.']);
end

if ~isequal(size(Kid),[2 2]) || numel(Tid) ~= 4
    error('Kid or Tid has unexpected dimensions.');
end

% Extract identified gains and time constants (+10% step)
K11 = Kid(1,1);
K12 = Kid(1,2);
K21 = Kid(2,1);
K22 = Kid(2,2);

T1 = Tid(1);
T2 = Tid(2);
T3 = Tid(3);
T4 = Tid(4);

fprintf('\n============================================\n');
fprintf('PROBLEM 4.6 - REALIZATION AND DISCRETIZATION\n');
fprintf('Identified model: +10%% step\n');
fprintf('Sample time: %.1f seconds\n',Ts);
fprintf('============================================\n');

fprintf('\nIdentified gain matrix:\n');
disp(Kid);

fprintf('Identified time constants [T1; T2; T3; T4]:\n');
disp(Tid);

%% Define the four transfer functions

transferFunctions(1).name = 'G11';
transferFunctions(1).num  = K11;
transferFunctions(1).den  = [T1 1];

transferFunctions(2).name = 'G12';
transferFunctions(2).num  = K12;
transferFunctions(2).den  = conv([T1 1],[T3 1]);

transferFunctions(3).name = 'G21';
transferFunctions(3).num  = K21;
transferFunctions(3).den  = conv([T2 1],[T4 1]);

transferFunctions(4).name = 'G22';
transferFunctions(4).num  = K22;
transferFunctions(4).den  = [T2 1];

%% Continuous-time realization and discretization

systems = struct();

for i = 1:4

    % Continuous-time observer canonical realization
    [Ac,Bc,Cc,Dc] = observerCanonical( ...
        transferFunctions(i).num, ...
        transferFunctions(i).den);

    % Exact zero-order-hold discretization
    [Ad,Bd,Cd,Dd] = discretizeZOH( ...
        Ac,Bc,Cc,Dc,Ts);

    % Store continuous-time matrices
    model.Ac = Ac;
    model.Bc = Bc;
    model.Cc = Cc;
    model.Dc = Dc;

    % Store discrete-time matrices
    model.Ad = Ad;
    model.Bd = Bd;
    model.Cd = Cd;
    model.Dd = Dd;

    systems.(transferFunctions(i).name) = model;

    %% Print discretized matrices

    fprintf('\n============================================\n');
    fprintf('%s - DISCRETE-TIME REALIZATION\n', ...
        transferFunctions(i).name);
    fprintf('============================================\n');

    disp('A_d ='); disp(Ad);
    disp('B_d ='); disp(Bd);
    disp('C_d ='); disp(Cd);
    disp('D_d ='); disp(Dd);

end

%% Store identified parameters

systems.K = Kid;
systems.T = Tid;
systems.Ts = Ts;

fprintf('\n============================================\n');
fprintf('PROBLEM 4.6 COMPLETED\n');
fprintf('============================================\n');

%% Observer canonical realization

function [A,B,C,D] = observerCanonical(num,den)

    % Normalize denominator to monic polynomial
    scale = den(1);

    den = den/scale;
    num = num/scale;

    % System order
    n = length(den)-1;

    % Pad numerator to n+1 coefficients
    num = [zeros(1,n+1-length(num)), num];

    % Extract coefficients
    a = den(2:end);
    b = num;

    % Observer canonical form
    A = zeros(n);

    A(:,1) = -a(:);

    for k = 1:n-1
        A(k,k+1) = 1;
    end

    B = (b(2:end) - a*b(1)).';

    C = [1 zeros(1,n-1)];

    D = b(1);

end

%% Zero-order-hold discretization

function [Ad,Bd,Cd,Dd] = discretizeZOH(Ac,Bc,Cc,Dc,Ts)

    n = size(Ac,1);
    m = size(Bc,2);

    % Augmented matrix
    M = [Ac, Bc;
         zeros(m,n), zeros(m,m)];

    % Matrix exponential
    Md = expm(M*Ts);

    % Extract discrete-time matrices
    Ad = Md(1:n,1:n);
    Bd = Md(1:n,n+1:n+m);

    Cd = Cc;
    Dd = Dc;

end


systems.Ts = Ts;


%% Compute MIMO Markov parameters H0, H1 and H2

% Mapping of the four SISO systems into the MIMO system
names = {'G11', 'G12';
         'G21', 'G22'};

% Initialize MIMO Markov parameter matrices
H = zeros(2,2,3);

for i = 1:2
    for j = 1:2

        % Retrieve corresponding discrete SISO realization
        sys = systems.(names{i,j});

        % H0 = D
        H(i,j,1) = sys.Dd;

        % H1 = C*B
        H(i,j,2) = sys.Cd * sys.Bd;

        % H2 = C*A*B
        H(i,j,3) = sys.Cd * sys.Ad * sys.Bd;

    end
end

%% Print MIMO Markov parameters

fprintf('\n============================================\n');
fprintf('MIMO MARKOV PARAMETERS\n');
fprintf('============================================\n');

fprintf('\nH_0 =\n');
disp(H(:,:,1));

fprintf('\nH_1 =\n');
disp(H(:,:,2));

fprintf('\nH_2 =\n');
disp(H(:,:,3));

% Store Markov parameters
systems.H = H;