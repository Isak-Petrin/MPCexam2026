function xdot = QuadrupleTankProcess(~, x, u, p)

% x = [h1; h2; h3; h4]
% u = [F1; F2]
% p = [a1; a2; a3; a4; A1; A2; A3; A4; g; gamma1; gamma2]


    % States
    h1 = x(1);
    h2 = x(2);
    h3 = x(3);
    h4 = x(4);

    % Inputs
    F1 = u(1);
    F2 = u(2);

    % Parameters
    a1 = p(1);
    a2 = p(2);
    a3 = p(3);
    a4 = p(4);

    A1 = p(5);
    A2 = p(6);
    A3 = p(7);
    A4 = p(8);

    g = p(9);

    gamma1 = p(10);
    gamma2 = p(11);

    % Four-tank model
    dh1 = gamma1/A1 * F1 ...
        + a3/A1 * sqrt(2*g*h3) ...
        - a1/A1 * sqrt(2*g*h1);

    dh2 = gamma2/A2 * F2 ...
        + a4/A2 * sqrt(2*g*h4) ...
        - a2/A2 * sqrt(2*g*h2);

    dh3 = (1-gamma2)/A3 * F2 ...
        - a3/A3 * sqrt(2*g*h3);

    dh4 = (1-gamma1)/A4 * F1 ...
        - a4/A4 * sqrt(2*g*h4);

    xdot = [dh1; dh2; dh3; dh4];

end