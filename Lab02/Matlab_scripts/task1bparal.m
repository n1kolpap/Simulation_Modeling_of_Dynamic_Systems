clear; clc;

m = 1.315; b = 0.225; k = 0.725;
theta_true = [b/m; k/m; 1/m]; % θ1, θ2, θ3

dt = 0.001; t = 0:dt:20;
u = 2.5 * sin(t)'; 
u_fun = @(tt) interp1(t, u, tt, 'linear', 'extrap');

% Αρχικές συνθήκες: [x1; x2; x1_hat; x2_hat; θ1_hat; θ2_hat; θ3_hat]
x0 = [0; 0; 0; 0; 0.2; 0.2; 0.2];

gamma = [0.1; 1; 1];

[t_out, X] = ode45(@(t, x) parallel_observer(t, x, u_fun, gamma), t, x0);


x1 = X(:,1); x2 = X(:,2);
x1hat = X(:,3); x2hat = X(:,4);
thetahat = X(:,5:7);

% Εκτιμήσεις παραμέτρων
bhat = thetahat(:,1) ./ thetahat(:,3);
khat = thetahat(:,2) ./ thetahat(:,3);
mhat = 1 ./ thetahat(:,3);

figure;
plot(t_out, mhat, 'b', t_out, bhat, 'r', t_out, khat, 'g'); hold on;
yline(m, '--b'); yline(b, '--r'); yline(k, '--g');
legend('m̂', 'b̂', 'k̂', 'm', 'b', 'k');
title('Εκτίμηση Παραμέτρων');
grid on;

figure;
plot(t_out, x1, 'b', t_out, x1hat, 'r');
legend('x', 'x̂'); 
grid on;

figure;
plot(t_out, x1 - x1hat);
title('Σφάλμα x - x̂');
grid on;

function dx = parallel_observer(t, x, u_fun, gamma)
    x1 = x(1); x2 = x(2);
    x1h = x(3); x2h = x(4);
    th1 = x(5); th2 = x(6); th3 = x(7);
    u = u_fun(t);

    % Σφάλματα
    e1 = x1 - x1h;
    e2 = x2 - x2h;


    dx1 = x2;
    dx2 = -0.225/1.315 * x2 - 0.725/1.315 * x1 + (1/1.315)*u;


    dx1h = x2h;
    dx2h = -th1 * x2h - th2 * x1h + th3 * u;


    dth1 = -gamma(1) * e2 * x2h;
    dth2 = -gamma(2) * e2 * x1h;
    dth3 =  gamma(3) * e2 * u;

    dx = [dx1; dx2; dx1h; dx2h; dth1; dth2; dth3];
end