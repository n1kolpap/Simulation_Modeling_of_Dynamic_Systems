clear; clc;

% παραμετροι του ελεγκτη
phi_inf = 0.001;
phi_0 = 50;
lambda = 5;
rho = 500;
k1 = 500;
k2 = 1000;

t = 0:0.001:20;

% επιθυμητη τροχια r_d(t) = (π/10) * sin(πt/20)
r_d_vec = (pi/10) * sin(pi * t / 20);
rd = @(tt) (pi/10) * sin(pi * tt / 20);


[t, x] = ode15s(@(t,x) odefun(t, x, phi_inf, phi_0, lambda, rho, k1, k2, rd), t, [0; 0]);

x1 = x(:,1); % r(t)
x2 = x(:,2); % ṙ(t)
phi = (phi_0 - phi_inf) * exp(-lambda * t) + phi_inf;


z1 = (x1 - r_d_vec') ./ phi;
a = -k1 * log((1 + z1) ./ (1 - z1));


figure();
plot(t, x1, 'b', 'LineWidth', 1.5); hold on;
plot(t, r_d_vec, 'r--', 'LineWidth', 1.2);
title('Αποκριση του συστηματος κλειστου βροχου');
xlabel('t [s]'); ylabel('r(t) [rad]');
legend('r(t)', 'r_d(t)');
grid on;

%  |r(t) - r_d(t)| < φ(t)
error = x1 - r_d_vec';
figure();
plot(t, error, 'k', 'LineWidth', 1.5); hold on;
plot(t, -phi, 'r--'); plot(t, phi, 'r--');
title('Έλεγχος συνθήκης |r(t) - r_d(t)| < φ(t)');
xlabel('t [s]'); ylabel('Σφάλμα r(t) - r_d(t)');
legend('Σφάλμα', '-φ(t)', 'φ(t)');
grid on;

%  |ṙ(t) - α(t)| < ρ
error2 = x2 - a;
rho_t = rho * ones(size(t));
figure();
plot(t, error2, 'b', 'LineWidth', 1.5); hold on;
plot(t, -rho_t, 'g--'); plot(t, rho_t, 'g--');
title('Έλεγχος συνθήκης |ṙ(t) - α(t)| < ρ');
xlabel('t [s]'); ylabel('Σφάλμα ṙ(t) - α(t)');
legend('Σφάλμα', '-ρ', 'ρ');
grid on;

function dx = odefun(t, x, phi_inf, phi_0, lambda, rho, k1, k2, rd)
    
    a1 = 1.315;
    a2 = 0.725;
    a3 = 0.225;
    b = 1.175;

    % Υπολογισμός ελεγκτή
    phi = (phi_0 - phi_inf) * exp(-lambda * t) + phi_inf;
    z1 = (x(1) - rd(t)) / phi;
    T1 = log((1 + z1) / (1 - z1));
    alpha = -k1 * T1;

    z2 = (x(2) - alpha) / rho;
    T2 = log((1 + z2) / (1 - z2));
    u = -k2 * T2;

    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = -a1*x(2) - a2*sin(x(1)) + a3*x(2)^2*sin(2*x(1)) + b*u;
end