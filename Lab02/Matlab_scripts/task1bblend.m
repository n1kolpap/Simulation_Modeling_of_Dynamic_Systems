clear; clc;

m = 1.315; b = 0.225; k = 0.725;
theta_true = [b/m; k/m; 1/m]; % θ1, θ2, θ3

dt = 0.001;
t = 0:dt:20;
u = 2.5 * sin(t)';
u_fun = @(tt) interp1(t, u, tt, 'linear', 'extrap');

% Αρχικές συνθήκες: [x1; x2; x1hat; x2hat; θ1hat; θ2hat; θ3hat]
x0 = [0; 0; 0; 0; 0.2; 0.2; 0.2];


gamma = [0.1; 0.5; 1];


[t_out, X] = ode45(@(t, x) mixed_observer(t, x, u_fun, gamma), t, x0);


x1     = X(:,1); 
x2     = X(:,2);
x1hat  = X(:,3); 
x2hat  = X(:,4);
thetahat = X(:,5:7);

% Εκτίμηση παραμέτρων
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
legend('x(t)', 'x̂(t)'); 


figure;
plot(t_out, x1 - x1hat);
title('Σφάλμα x - x̂');
xlabel('Χρόνος'); ylabel('Σφάλμα'); grid on;


function dx = mixed_observer(t, x, u_fun, gamma)
   
    x1 = x(1); x2 = x(2);
    x1hat = x(3); x2hat = x(4);
    th1 = x(5); th2 = x(6); th3 = x(7);
    u = u_fun(t);

    % Σφάλματα 
    e1 = x1 - x1hat;
    e2 = x2 - x2hat;

 
    dx1 = x2;
    dx2 = -0.225/1.315 * x2 - 0.725/1.315 * x1 + (1/1.315)*u;

    % Εκτιμητής
    lambda1 = 10; %θm1
    lambda2 = 10; %θm2
    dx1hat = x2hat + lambda1 * e1;
    dx2hat = -th1 * x2hat - th2 * x1hat + th3 * u + lambda2 * e2;

    dth1 = -gamma(1) * e2 * x2hat;
    dth2 = -gamma(2) * e2 * x1hat;
    dth3 =  gamma(3) * e2 * u;

    dx = [dx1; dx2; dx1hat; dx2hat; dth1; dth2; dth3];
end
