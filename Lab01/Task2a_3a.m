clear;
clc;
close all;

m = 0.75;
L = 1.25;
c = 0.15;
g = 9.81;

x0 = [0; 0];                    
t_full = 0:0.001:20;            
u = @(t) 4 * sin(2 * t);        

f = @(t, x) [x(2); (1/(m*L^2)) * (u(t) - c*x(2) - m*g*L*x(1))];

[~, x] = ode45(f, t_full, x0);       
y_full = x(:,1);                     
dy_full = x(:,2);                    

% Δειγματοληψια
Ts = 0.1;
t = 0:Ts:20;
y = interp1(t_full, y_full, t)';
dy = interp1(t_full, dy_full, t)';
u_sampled = 4 * sin(2 * t)';

% Φιλτρο 
s = tf('s');
sys1 = s / (s + 1);          % s*q / (s+1)
sys2 = 1 / (s + 1);          % q / (s+1)
sys3 = 1 / (s + 1);          % u / (s+1)


z1 = lsim(sys1, y, t);       % s*q / (s+1)
z2 = lsim(sys2, y, t);       % q / (s+1)
z3 = lsim(sys3, u_sampled, t);  % u / (s+1)

% Least Squares
Z = [z1 z2 z3];   % Φ
Qdot = dy;        % q̇(t) 
a = Z' * Z;
b = Qdot' * Z;
theta = b / a;

% Εκτίμηση παραμέτρων
ML2 = 1 / theta(3);                         % mL^2
chat = (1 - theta(1)) * ML2;                % ĉ
mgL = -theta(2) * ML2;                      % mgL
mhat = mgL / (g * L);                       % m̂
lhat = sqrt(ML2 / mhat);                    % L̂

fprintf('\nΕκτιμήσεις παραμέτρων\n');
fprintf('m̂ = %.4f kg\n', mhat);
fprintf('L̂ = %.4f m\n', lhat);
fprintf('ĉ = %.4f N·m·s\n', chat);

fhat = @(t, x) [x(2); (1/(mhat*lhat^2)) * (u(t) - chat*x(2) - mhat*g*lhat*x(1))];
[~, xhat] = ode45(fhat, t_full, x0);
yhat = xhat(:,1);   % q̂(t)

figure;
plot(t_full, y_full, 'b', t_full, yhat, 'r--');
legend('q(t)', 'q̂(t)');
xlabel('t [sec]');
ylabel('q(t)');
grid on;

% calculate error and plot(t, e)
e = y_full - yhat;
figure;
plot(t_full, e)
title('Σφάλμα e_q(t) = q(t) - q̂(t)');
xlabel('t [s]');
ylabel('e_q [m]');
grid on;
% ΘΕΜΑ 3α
noise_level = 0.1;

y_noisy = y + noise_level * randn(size(y));            % q(t) με θόρυβο
dy_noisy = dy + noise_level * randn(size(dy));         % q̇(t) με θόρυβο
u_noisy = u_sampled + noise_level * randn(size(u_sampled));  % u(t) με θόρυβο

z1_noisy = lsim(sys1, y_noisy, t);  
z2_noisy = lsim(sys2, y_noisy, t);  
z3_noisy = lsim(sys3, u_noisy, t);  

Z_noisy = [z1_noisy z2_noisy z3_noisy];  % νέο Φ
Qdot_noisy = dy_noisy;                   % νέο q̇(t)

a_noisy = Z_noisy' * Z_noisy;
b_noisy = Qdot_noisy' * Z_noisy;
theta_noisy = b_noisy / a_noisy;

ML2_n = 1 / theta_noisy(3);
chat_n = (1 - theta_noisy(1)) * ML2_n;
mgL_n = -theta_noisy(2) * ML2_n;
mhat_n = mgL_n / (g * L);
lhat_n = sqrt(ML2_n / mhat_n);

% Εκτύπωση συγκρίσεων
fprintf('\nΕκτιμήσεις χωρίς θόρυβο\n');
fprintf('m̂ = %.4f kg\n', mhat);
fprintf('L̂ = %.4f m\n', lhat);
fprintf('ĉ = %.4f N·m·s\n', chat);

fprintf('\nΕκτιμήσεις με θόρυβο\n');
fprintf('m̂ (noisy) = %.4f kg\n', mhat_n);
fprintf('L̂ (noisy) = %.4f m\n', lhat_n);
fprintf('ĉ (noisy) = %.4f N·m·s\n', chat_n);

fprintf('\nΑπόκλιση\n');
fprintf('Δm = %.4f\n', abs(mhat - mhat_n));
fprintf('ΔL = %.4f\n', abs(lhat - lhat_n));
fprintf('Δc = %.4f\n', abs(chat - chat_n));

yhat_noisy = Z_noisy * theta_noisy';

figure;
plot(t, y, 'b', t, yhat_noisy, 'r--');
legend('q(t)', 'q̂(t) with noise');
title('Σύγκριση q(t) με θόρυβο');
xlabel('t [s]');
ylabel('q(t)');
grid on;

% Σφάλμα εξόδου με θόρυβο
e_noise = y - yhat_noisy;
figure;
plot(t, e_noise);
title('Σφάλμα με θόρυβο: e_q(t) = q(t) - q̂(t)');
xlabel('t [s]');
ylabel('e_q(t)');
grid on;