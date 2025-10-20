clear;
clc;
close all;

% Πραγματικές τιμές
m = 0.75;
L = 1.25;
c = 0.15;
g = 9.81;
Ts = 0.1;

x0 = [0; 0];
t_full = 0:0.001:20;
f_real = @(t, x, u) [x(2); (1/(m*L^2)) * (u - c*x(2) - m*g*L*x(1))];

% Πλάτη εισόδου A0 προς μελέτη
A_values = [0.5 1 2 4 6 8 10];
errors_m = zeros(size(A_values));
errors_L = zeros(size(A_values));
errors_c = zeros(size(A_values));

for i = 1:length(A_values)
    A = A_values(i);
    u_fun = @(t) A * sin(2 * t);

    f = @(t, x) f_real(t, x, u_fun(t));
    [~, x_full] = ode45(f, t_full, x0);
    y_full = x_full(:,1);
    dy_full = x_full(:,2);

    t = 0:Ts:20;
    y = interp1(t_full, y_full, t)';
    dy = interp1(t_full, dy_full, t)';
    u = A * sin(2 * t)';

    % Φίλτρα Λ(s) = s + 1
    s = tf('s');
    sys1 = s / (s + 1);
    sys2 = 1 / (s + 1);
    sys3 = 1 / (s + 1);

    z1 = lsim(sys1, y, t);  % s*q / (s+1)
    z2 = lsim(sys2, y, t);  % q / (s+1)
    z3 = lsim(sys3, u, t);  % u / (s+1)

    Z = [z1 z2 z3];
    Qdot = dy;

    theta = (Qdot' * Z) / (Z' * Z);

    % Εκτίμηση παραμέτρων
    ML2 = 1 / theta(3);
    chat = (1 - theta(1)) * ML2;
    mgL = -theta(2) * ML2;
    mhat = mgL / (g * L);
    lhat = sqrt(ML2 / mhat);

    % Σφάλματα εκτίμησης
    errors_m(i) = abs(m - mhat);
    errors_L(i) = abs(L - lhat);
    errors_c(i) = abs(c - chat);
end

figure;
plot(A_values, errors_m, '-o', 'DisplayName','σφάλμα m');
hold on;
plot(A_values, errors_L, '-s', 'DisplayName','σφάλμα L');
plot(A_values, errors_c, '-^', 'DisplayName','σφάλμα c');
xlabel('Πλάτος εισόδου A₀');
ylabel('Σφάλμα εκτίμησης');
legend;
grid on;
