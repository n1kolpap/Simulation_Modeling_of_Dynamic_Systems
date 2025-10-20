clear; clc;

% Παράμετροι ελέγχου
phi_inf = 0.01;
phi_0 = 100;
lambda = 8;
rho = 100;
k1 = 10;
k2 = 20;

t = 0:0.001:20;
r1 = (pi/10) * sin(pi * t / 20);
rd = @(t) (pi/10) * sin(pi * t / 20);


[t,x] = ode15s(@(t,x) odefun(t,x,phi_inf,phi_0,lambda,rho,k1,k2,rd), t, [0;0]);

phi = (phi_0 - phi_inf) * exp(-lambda * t) + phi_inf;
z1 = (x(:,1) - r1') ./ phi;
T1 = log((1 + z1) ./ (1 - z1));
a = -k1 * T1;
z2 = (x(:,2) - a) ./ rho;
T2 = log((1 + z2) ./ (1 - z2));
u = -k2 * T2;


t_data = t;
x2 = @(t) interp1(t_data, x(:,2), t, 'linear', 'extrap');
x1 = @(t) interp1(t_data, x(:,1), t, 'linear', 'extrap');
u  = @(t) interp1(t_data, u, t, 'linear', 'extrap');


gamma = [1800,900,1400,600];
[t,z] = ode15s(@(t,Z) odefun2(t,Z,x1,x2,u,gamma), t, [0;0;0.01;0.01;0.25;0.01]);

x_hat = z(:,1);
a1_est = z(:,3);
a2_est = z(:,4);
a3_est = z(:,5);
b_est  = z(:,6);
error  = x(:,1) - x_hat;

a1 = 1.315 * ones(length(t),1);
a2 = 0.725 * ones(length(t),1);
a3 = 0.225 * ones(length(t),1);
b  = 1.175 * ones(length(t),1);

figure();
plot(t, x(:,1), 'b', t, x_hat, 'r');
legend('r(t)', 'r_{est}(t)');
xlabel('t [s]'); ylabel('Γωνία [rad]');
title('Πραγματική και εκτιμημένη γωνία r(t) και r_{est}(t)');
grid on;


figure();
plot(t, error, 'k');
xlabel('t [s]'); ylabel('Σφάλμα [rad]');
title('Σφάλμα εκτίμησης: r(t) - r_{est}(t)');
grid on;


figure();
plot(t, a1_est, 'b', t, a1, '--k');
legend('a1 est','a1 true');
xlabel('t [s]');
title('Εκτίμηση παραμέτρου a1');
grid on;


figure();
plot(t, a2_est, 'r', t, a2, '--k');
legend('a2 est','a2 true');
xlabel('t [s]');
title('Εκτίμηση παραμέτρου a2');
grid on;


figure();
plot(t, a3_est, 'g', t, a3, '--k');
legend('a3 est','a3 true');
xlabel('t [s]');
title('Εκτίμηση παραμέτρου a3');
grid on;

figure();
plot(t, b_est, 'm', t, b, '--k');
legend('b est','b true');
xlabel('t [s]');
title('Εκτίμηση παραμέτρου b');
grid on;

function dx = odefun(t,x,phi_inf,phi_0,lambda,rho,k1,k2,rd)
    a1 = 1.315; a2 = 0.725; a3 = 0.225; b = 1.175;
    d = 0.15 * sin(0.5 * t); 
    phi = (phi_0 - phi_inf) * exp(-lambda * t) + phi_inf;
    z1 = (x(1) - rd(t)) / phi;
    T1 = log((1 + z1) / (1 - z1));
    alpha = -k1 * T1;
    z2 = (x(2) - alpha) / rho;
    T2 = log((1 + z2) / (1 - z2));
    u = -k2 * T2;
    dx = [x(2); -a1*x(2) - a2*sin(x(1)) + a3*x(2)^2*sin(2*x(1)) + b*u + d]; 
end

function dZ = odefun2(t,Z,x1,x2,u,gamma)
    x = Z(1:2); theta = Z(3:6);
    dx1 = x2(t) + 100*(x1(t) - x(1));
    dx2 = -theta(1)*x2(t) - theta(2)*sin(x1(t)) + theta(3)*x2(t)^2*sin(2*x1(t)) + theta(4)*u(t) + 100*(x2(t) - x(2));
    dth1 = -gamma(1)*(x2(t)-x(2))*x2(t);
    dth2 = -gamma(2)*sin(x1(t))*(x2(t)-x(2));
    dth3 =  gamma(3)*x2(t)^2*sin(2*x1(t))*(x2(t)-x(2));
    dth4 =  gamma(4)*u(t)*(x2(t)-x(2));
    dZ = [dx1; dx2; dth1; dth2; dth3; dth4];
end
