clear; clc;

m_real = 1.315;
b_real = 0.225;
k_real = 0.725;


dt = 0.001;
t = 0:dt:20;


u_type = 2;
%u_type = 1;
if u_type == 1
    u = 2.5 * ones(length(t), 1);
else
    u = 2.5 * sin(t)';
end
u_fun = @(tt) interp1(t, u, tt);


x0 = [0; 0];
[t_out, x_out] = ode45(@(t, x) real_sys(t, x, u_fun, m_real, b_real, k_real), t, x0);
x = x_out(:,1);
xdot = x_out(:,2);

% φιλτρα
a_m = 1;
phi1 = lsim(tf([-1 0], [1 a_m]), xdot, t);
phi2 = lsim(tf([0 -1], [1 a_m]), x, t);
phi3 = lsim(tf([0 1], [1 a_m]), u, t);

phi1_fun = @(tt) interp1(t, phi1, tt);
phi2_fun = @(tt) interp1(t, phi2, tt);
phi3_fun = @(tt) interp1(t, phi3, tt);
y_dot_fun = @(tt) interp1(t, xdot, tt);

% εκτιμηση παραμετρων
gamma = diag([10, 10, 10]);
theta0 = [-0.1; 0.1; 0.1];
[t_theta, thetahat] = ode45(@(t, theta) gradient_update(t, theta, phi1_fun, phi2_fun, phi3_fun, y_dot_fun, gamma), t, theta0);
thetahat = thetahat';

a1_hat = thetahat(1,:) + a_m;
a2_hat = thetahat(2,:);
b1_hat = thetahat(3,:);

m_hat = 1 ./ b1_hat;
b_hat = a1_hat .* m_hat;
k_hat = a2_hat .* m_hat;

% πραγματικες τιμες
m_true = m_real * ones(size(t));
b_true = b_real * ones(size(t));
k_true = k_real * ones(size(t));

% γραφικες παραστασεις
figure;
plot(t, m_hat, 'b'); hold on;
plot(t, b_hat, 'r');
plot(t, k_hat, 'g');
plot(t, m_true, '--b'); 
plot(t, b_true, '--r'); 
plot(t, k_true, '--g');
legend('m_hat.', 'b_hat.', 'k_hat.', 'm', 'b', 'k');
title('Εκτίμηση Παραμέτρων');
grid on;

% x̂(t)
phi_mat = [phi1'; phi2'; phi3'];
xhat_dot = sum(thetahat .* phi_mat);
xhat = cumtrapz(t, xhat_dot);

% σφαλμα μοντελοποιησης
error = x - xhat';

figure;
plot(t, x, 'b'); hold on;
plot(t, xhat, 'r');
legend('x', 'x̂');
grid on;

figure;
plot(t, error);
title('Σφάλμα x - x̂');
grid on;

% συναρτησεις
function dx = real_sys(t, x, u_fun, m, b, k)
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = -(b/m)*x(2) - (k/m)*x(1) + (1/m)*u_fun(t);
end

function dtheta = gradient_update(t, theta, phi1, phi2, phi3, ydot, gamma)
    phi_vec = [phi1(t); phi2(t); phi3(t)];
    e = ydot(t) - theta' * phi_vec;
    dtheta = gamma * e * phi_vec;
end