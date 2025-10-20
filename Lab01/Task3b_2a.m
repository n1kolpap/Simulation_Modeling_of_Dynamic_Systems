clear;
clc;

% Πραγματικές τιμές
m = 0.75;
L = 1.25;
c = 0.15;
g = 9.81;

x0 = [0; 0];
t_full = 0:0.001:20;
u_fun = @(t) 4 * sin(2 * t);
f = @(t, x) [x(2); (1/(m*L^2)) * (u_fun(t) - c*x(2) - m*g*L*x(1))];
[~, x_full] = ode45(f, t_full, x0);
y_full = x_full(:,1);
dy_full = x_full(:,2);

% Ts τιμές 
Ts_values = [0.01 0.02 0.05 0.1 0.2 0.5];
errors_m = zeros(size(Ts_values));
errors_L = zeros(size(Ts_values));
errors_c = zeros(size(Ts_values));

ML2_real = m * L^2;
chat_real = c;
mgL_real = m * g * L;

for i = 1:length(Ts_values)
    Ts = Ts_values(i);
    t = 0:Ts:20;
    
    y = interp1(t_full, y_full, t)';
    dy = interp1(t_full, dy_full, t)';
    u = 4 * sin(2 * t)';

    s = tf('s');
    sys1 = s / (s + 1);
    sys2 = 1 / (s + 1);
    sys3 = 1 / (s + 1);

    z1 = lsim(sys1, y, t);
    z2 = lsim(sys2, y, t);
    z3 = lsim(sys3, u, t);

    Z = [z1 z2 z3];
    Qdot = dy;

    theta = (Qdot' * Z) / (Z' * Z);

    ML2 = 1 / theta(3);
    chat = (1 - theta(1)) * ML2;
    mgL = -theta(2) * ML2;
    mhat = mgL / (g * L);
    lhat = sqrt(ML2 / mhat);

    errors_m(i) = abs(m - mhat);
    errors_L(i) = abs(L - lhat);
    errors_c(i) = abs(c - chat);
end

figure;
plot(Ts_values, errors_m, '-o', 'DisplayName','σφάλμα m');
hold on;
plot(Ts_values, errors_L, '-s', 'DisplayName','σφάλμα L');
plot(Ts_values, errors_c, '-^', 'DisplayName','σφάλμα c');
xlabel('Ts [sec]');
ylabel('Σφάλμα Εκτίμησης');
legend;
grid on;

