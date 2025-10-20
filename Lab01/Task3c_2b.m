clear;
clc;
close all;

m = 0.75;
L = 1.25;
c = 0.15;
g = 9.81;
Ts = 0.1;

x0 = [0; 0];
t_full = 0:0.001:20;
f_real = @(t, x, u) [x(2); (1/(m*L^2)) * (u - c*x(2) - m*g*L*x(1))];

A_values = [0.5 1 2 4 6 8 10];
errors_m = zeros(size(A_values));
errors_L = zeros(size(A_values));
errors_c = zeros(size(A_values));

ML2_real = m * L^2;
chat_real = c;
mgL_real = m * g * L;

for i = 1:length(A_values)
    A = A_values(i);
    u = @(t) A * sin(2 * t);

    f = @(t, x) f_real(t, x, u(t));
    [~, x] = ode45(f, t_full, x0);
    y_full = x(:,1);

    t = 0:Ts:20;
    y = interp1(t_full, y_full, t)';
    u_sampled = A * sin(2 * t)';

    % Φίλτρα (Λ(s) = s^2 + 2s + 1)
    l1 = 2; l2 = 1;
    s = tf('s');
    sys1 = tf([-1 0], [1 l1 l2]);   % -s*q / Λ(s)
    sys2 = tf(-1, [1 l1 l2]);       % -q / Λ(s)
    sys3 = tf(1, [1 l1 l2]);        % u / Λ(s)

    z1 = lsim(sys1, y, t);
    z2 = lsim(sys2, y, t);
    z3 = lsim(sys3, u_sampled, t);

    Z = [z1 z2 z3];  % Φ
    y_interp = y;    % q(t)
    theta = (y_interp' * Z) / (Z' * Z);  % θ = y^T * Z / (Z^T * Z)

    % Εκτίμηση παραμέτρων
    ML2 = 1 / theta(3);
    chat = (theta(1) + l1) * ML2;
    mgL = (theta(2) + l2) * ML2;
    mhat = mgL / (g * L);
    lhat = sqrt(ML2 / mhat);

    % Υπολογισμός σφαλμάτων
    errors_m(i) = abs(m - mhat);
    errors_L(i) = abs(L - lhat);
    errors_c(i) = abs(c - chat);
end

% Γραφήματα
figure;
plot(A_values, errors_m, '-o', 'DisplayName','σφάλμα m');
hold on;
plot(A_values, errors_L, '-s', 'DisplayName','σφάλμα L');
plot(A_values, errors_c, '-^', 'DisplayName','σφάλμα c');
xlabel('Πλάτος εισόδου A₀');
ylabel('Σφάλμα εκτίμησης');
title('Εκτίμηση Παραμέτρων vs Πλάτος Εισόδου A₀');
legend;
grid on;
