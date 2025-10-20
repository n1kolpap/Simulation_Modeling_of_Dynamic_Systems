clear;
clc;
close all;

m = 0.75;
L = 1.25;
c = 0.15;
g = 9.81;

x0 = [0; 0];                    
t = 0:0.001:20;                 
u = @(t) 4 * sin(2 * t);        

f = @(t, x) [x(2); (1/(m*L^2)) * (u(t) - c*x(2) - m*g*L*x(1))];

[t, x] = ode45(f, t, x0);       
y = x(:,1);                     % q(t)
dy = x(:,2);                    % q̇(t)

% Δειγματοληψια
Ts = 0.1;
t_sampled = 0:Ts:20;
y = interp1(t, y, t_sampled)';
u_sampled = 4 * sin(2 * t_sampled)';

% Λ(s) = s^2 + 2s + 1
l1 = 2;
l2 = 1;
sys1 = tf([-1 0], [1 l1 l2]);
sys2 = tf(-1, [1 l1 l2]);
sys3 = tf(1, [1 l1 l2]);

[~, y_full] = ode45(f, t_sampled, x0);
y_interp = y_full(:,1);   
t = t_sampled;       

z1 = lsim(sys1, y_interp, t);    % -s*q / Λ(s)
z2 = lsim(sys2, y_interp, t);    % -q / Λ(s)
z3 = lsim(sys3, u_sampled, t);   % u / Λ(s)

% Least Squares
Z = [z1 z2 z3];   % Φ
a = Z' * Z;
b = y_interp' * Z;
theta = b / a;

% Εκτιμηση παραμέτρων
ML2 = 1 / theta(3);                         % mL^2
chat = (theta(1) + l1) * ML2;               % ĉ
mgL = (theta(2) + l2) * ML2;                % mgL
mhat = mgL / (g * L);                       % m̂
lhat = sqrt(ML2 / mhat);                    % L̂

fprintf('\nΕκτιμήσεις παραμέτρων\n');
fprintf('m̂ = %.4f kg\n', mhat);
fprintf('L̂ = %.4f m\n', lhat);
fprintf('ĉ = %.4f N·m·s\n', chat);

yhat = Z * theta';
figure;
plot(t, y_interp, 'b', t, yhat, 'r--');
legend('q(t)', 'q̂(t)');
xlabel('t [sec]');
ylabel('q(t)');
grid on;

e = y - yhat;
figure;
plot(t, e)
title('e(t) - t');
xlabel('t [s]');
ylabel('e [m]');
grid on;

% ΘΕΜΑ 3α
noise_level = 0.1;

y_noisy = y_interp + noise_level * randn(size(y_interp));  % θόρυβος στο q(t)
u_noisy = u_sampled + noise_level * randn(size(u_sampled)); 

z1_noisy = lsim(sys1, y_noisy, t);  
z2_noisy = lsim(sys2, y_noisy, t);  
z3_noisy = lsim(sys3, u_noisy, t);  

Z_noisy = [z1_noisy z2_noisy z3_noisy];
a_noisy = Z_noisy' * Z_noisy;
b_noisy = y_noisy' * Z_noisy;
theta_noisy = b_noisy / a_noisy;

% Εκτίμηση παραμέτρων με θόρυβο
ML2_n = 1 / theta_noisy(3);
chat_n = (theta_noisy(1) + l1) * ML2_n;
mgL_n = (theta_noisy(2) + l2) * ML2_n;
mhat_n = mgL_n / (g * L);
lhat_n = sqrt(ML2_n / mhat_n);

fprintf('\nΕκτιμήσεις χωρίς θόρυβο\n');
fprintf('m̂ = %.4f kg\n', mhat);
fprintf('L̂ = %.4f m\n', lhat);
fprintf('ĉ = %.4f N·m·s\n', chat);

fprintf('\nΕκτιμήσεις με θόρυβο\n');
fprintf('m̂ (noisy) = %.4f kg\n', mhat_n);
fprintf('L̂ (noisy) = %.4f m\n', lhat_n);
fprintf('ĉ (noisy) = %.4f N·m·s\n', chat_n);

% Σφάλματα
fprintf('\nΑπόκλιση\n');
fprintf('Δm = %.4f\n', abs(mhat - mhat_n));
fprintf('ΔL = %.4f\n', abs(lhat - lhat_n));
fprintf('Δc = %.4f\n', abs(chat - chat_n));

yhat_noisy = Z_noisy * theta_noisy';

figure;
plot(t, y_interp, 'b', t, yhat_noisy, 'r--');
legend('q(t)', 'q̂(t) with noise');
title('Σύγκριση q(t) με θόρυβο');
xlabel('t [s]');
ylabel('q(t)');
grid on;

% Σφάλμα με θόρυβο
e_noise = y_interp - yhat_noisy;
figure;
plot(t, e_noise);
title('Σφάλμα με θόρυβο: e_q(t) = q(t) - q̂(t)');
xlabel('t [s]');
ylabel('e_q(t)');
grid on;