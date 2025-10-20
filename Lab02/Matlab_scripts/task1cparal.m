clear; clc;

% Πραγματικές παράμετροι
m = 1.315; b = 0.225; k = 0.725;

% Χρονικά δεδομένα
dt = 0.0001;
t = 0:dt:20;
u = 2.5 * sin(t)';
u_fun = @(tt) interp1(t, u, tt, 'linear', 'extrap');

% Πλάτος και συχνότητα θορύβου
eta0 = 0.25; 
f0 = 20;
noise = eta0 * sin(2*pi*f0*t); % η(t)

% Αρχικές συνθήκες
x0 = [0; 0]; % Πραγματικό σύστημα
[t_real, x_real] = ode15s(@(t,x) true_system(t, x, u_fun), t, x0);

x1 = x_real(:,1); % καθαρό x(t)
x2 = x_real(:,2);

x1_noisy = x1 + noise'; % μετρημένο x(t) με θόρυβο
x1_noisy_fun = @(tt) interp1(t, x1_noisy, tt, 'linear', 'extrap');
x2_fun = @(tt) interp1(t, x2, tt, 'linear', 'extrap');

gamma = [0.1; 1; 1];
xhat0 = [0; 0; 0.2; 0.2; 0.2]; % [x1_hat; x2_hat; θ1; θ2; θ3]

[t_out, Z] = ode15s(@(t,z) parallel_observer_mixed(t,z,u_fun,x1_noisy_fun,x2_fun,gamma), t, xhat0);

x1_hat = Z(:,1);
theta1 = Z(:,3); theta2 = Z(:,4); theta3 = Z(:,5);

m_hat = 1 ./ theta3;
b_hat = theta1 ./ theta3;
k_hat = theta2 ./ theta3;

figure;
plot(t_out, m_hat, 'b', t_out, b_hat, 'r', t_out, k_hat, 'g'); hold on;
yline(m, '--b'); yline(b, '--r'); yline(k, '--g');
legend('m̂','b̂','k̂','m','b','k');
title('Εκτίμηση παραμέτρων');
grid on;

figure;
plot(t, x1, 'b', t_out, x1_hat, 'r');
legend('x(t)','x̂(t)');
grid on;

figure;
plot(t_out, x1_noisy(1:length(t_out)) - x1_hat);
title('Σφάλμα x(t) - x̂(t) με θόρυβο');
grid on;


function dx = true_system(t,x,u_fun)
    m = 1.315; b = 0.225; k = 0.725;
    u = u_fun(t);
    dx = [x(2);
         -(b/m)*x(2) - (k/m)*x(1) + (1/m)*u];
end

function dz = parallel_observer_mixed(t,z,u_fun,x1_meas_fun,x2_fun,gamma)
    x1_hat = z(1); x2_hat = z(2);
    th1 = z(3); th2 = z(4); th3 = z(5);

    u = u_fun(t);
    x1_meas = x1_meas_fun(t); % με θόρυβο
    x2 = x2_fun(t); % καθαρό

    e1 = x1_meas - x1_hat;
    e2 = x2 - x2_hat;

    dx1h = x2_hat;
    dx2h = -th1 * x2_hat - th2 * x1_hat + th3 * u;

    dth1 = -gamma(1) * e2 * x2_hat;
    dth2 = -gamma(2) * e2 * x1_hat;
    dth3 =  gamma(3) * e2 * u;

    dz = [dx1h; dx2h; dth1; dth2; dth3];
end
