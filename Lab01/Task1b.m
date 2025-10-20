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

y = x(:,1);
dy = x(:,2);

figure(1);
plot(t, y);
xlabel('t [s]');
ylabel('q(t) [rad]');
title('q(t) - t');
grid on;

figure(2);
plot(t, dy);
xlabel('t [s]');
ylabel('dq/dt(t) [rad/s]');
title('dq/dt(t) - t');
grid on;
