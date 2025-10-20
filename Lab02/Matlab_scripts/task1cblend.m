clear; clc;


m = 1.315; b = 0.225; k = 0.725;
dt = 0.001;
t = 0:dt:20;
u = 2.5 * sin(t)';
u_fun = @(tt) interp1(t, u, tt, 'linear', 'extrap');
gamma = [0.1; 0.5; 1];
lambda1 = 10; lambda2 = 10;


n0 = 1; f0 = 20;
noise = n0 * sin(2*pi*f0*t)';
[~, x_real] = ode15s(@(t,x) real_system(t,x,u_fun), t, [0; 0]);
x1 = x_real(:,1); x2 = x_real(:,2);
x1_noisy = x1 + noise;

x1_fun = @(tt) interp1(t, x1_noisy, tt, 'linear', 'extrap');
x2_fun = @(tt) interp1(t, x2, tt, 'linear', 'extrap');
x0 = [0; 0; 0.2; 0.2; 0.2];

[t_out, Z] = ode15s(@(t,z) mixed_noise_observer(t,z,u_fun,x1_fun,x2_fun,gamma,lambda1,lambda2), t, x0);
x1hat = Z(:,1); x2hat = Z(:,2);
theta1 = Z(:,3); theta2 = Z(:,4); theta3 = Z(:,5);

mhat = 1 ./ theta3;
bhat = theta1 ./ theta3;
khat = theta2 ./ theta3;
error = x1 - x1hat;

% Σφάλματα εκτίμησης παραμέτρων 
e_m = abs(mhat - m) / m;
e_b = abs(bhat - b) / b;
e_k = abs(khat - k) / k;

figure;
plot(t, e_m, 'b', t, e_b, 'r', t, e_k, 'g');
legend('e_m', 'e_b', 'e_k');
title('Σφάλμα Εκτίμησης Παραμέτρων με \eta_0 = 0.25');
xlabel('Χρόνος [s]');
ylabel('Σχετικό Σφάλμα');
grid on;


figure;
plot(t, error);
title('Σφάλμα e(t) = x(t) - x̂(t)');
grid on;


n0_vals = 0:0.1:1;
E_m = zeros(length(t), length(n0_vals));
E_b = zeros(length(t), length(n0_vals));
E_k = zeros(length(t), length(n0_vals));

for i = 1:length(n0_vals)
    n0 = n0_vals(i);
    noise = n0 * sin(2*pi*f0*t)';
    [~, x_real] = ode15s(@(t,x) real_system(t,x,u_fun), t, [0; 0]);
    x1 = x_real(:,1); x2 = x_real(:,2);
    x1_noisy = x1 + noise;

    x1_fun = @(tt) interp1(t, x1_noisy, tt, 'linear', 'extrap');
    x2_fun = @(tt) interp1(t, x2, tt, 'linear', 'extrap');
    x0 = [0; 0; 0.2; 0.2; 0.2];

    [~, Z] = ode15s(@(t,z) mixed_noise_observer(t,z,u_fun,x1_fun,x2_fun,gamma,lambda1,lambda2), t, x0);

    theta1 = Z(:,3); theta2 = Z(:,4); theta3 = Z(:,5);
    mhat = 1 ./ theta3;
    bhat = theta1 ./ theta3;
    khat = theta2 ./ theta3;

    E_m(:,i) = abs(mhat - m);
    E_b(:,i) = abs(bhat - b);
    E_k(:,i) = abs(khat - k);
end


figure;
plot(t, E_m); legend(arrayfun(@(x) sprintf('\\eta_0=%.1f',x), n0_vals, 'UniformOutput', false));
title('Σφάλμα m̂ για διαφορετικά \eta_0'); grid on;

figure;
plot(t, E_b); legend(arrayfun(@(x) sprintf('\\eta_0=%.1f',x), n0_vals, 'UniformOutput', false));
title('Σφάλμα b̂ για διαφορετικά \eta_0'); grid on;

figure;
plot(t, E_k); legend(arrayfun(@(x) sprintf('\\eta_0=%.1f',x), n0_vals, 'UniformOutput', false));
title('Σφάλμα k̂ για διαφορετικά \eta_0'); grid on;

function dx = real_system(t,x,u_fun)
    m = 1.315; b = 0.225; k = 0.725;
    u = u_fun(t);
    dx = [x(2); -(b/m)*x(2) - (k/m)*x(1) + (1/m)*u];
end

function dz = mixed_noise_observer(t,z,u_fun,x1_fun,x2_fun,gamma,lambda1,lambda2)
    x1hat = z(1); x2hat = z(2);
    th1 = z(3); th2 = z(4); th3 = z(5);
    u = u_fun(t);

    x1_meas = x1_fun(t);
    x2_meas = x2_fun(t);

    e1 = x1_meas - x1hat;
    e2 = x2_meas - x2hat;

    dx1h = x2hat + lambda1 * e1;
    dx2h = -th1 * x2hat - th2 * x1hat + th3 * u + lambda2 * e2;

    dth1 = -gamma(1) * e2 * x2hat;
    dth2 = -gamma(2) * e2 * x1hat;
    dth3 =  gamma(3) * e2 * u;

    dz = [dx1h; dx2h; dth1; dth2; dth3];
end