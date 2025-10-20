T = 20; dt = 0.01; N = T/dt;
t = linspace(0, T, N);


A = [-2.15  0.25;
          -0.75 -2.00];
B = [0; 1.5];


%omega_bar = 0.5;
omega_bar=5;


gamma1 = 5; gamma2 = 5; sigma = 0.1;
theta_m1 = 10; theta_m2 = 10;


x = zeros(2, N); x_hat = zeros(2, N);
theta1_hat = [-2; 0; 0];
theta2_hat = [-1; 0; 0];
theta1_hist = zeros(3, N); theta2_hist = zeros(3, N);


u = sin(2 * t);

for k = 1:N-1
    phi = [x(1,k); x(2,k); u(k)];

  
    omega = omega_bar * (2*rand(2,1)-1);


    dx = A * x(:,k) + B * u(k) + omega;
    x(:,k+1) = x(:,k) + dt * dx;

    
    x1_hat_dot = phi' * theta1_hat + theta_m1 * (x(1,k) - x_hat(1,k));
    x2_hat_dot = phi' * theta2_hat + theta_m2 * (x(2,k) - x_hat(2,k));
    x_hat(1,k+1) = x_hat(1,k) + dt * x1_hat_dot;
    x_hat(2,k+1) = x_hat(2,k) + dt * x2_hat_dot;

   
    ex1 = x(1,k) - x_hat(1,k);
    ex2 = x(2,k) - x_hat(2,k);

    % Εκτίμηση παραμέτρων με σ-τροποποίηση
    dtheta1 = gamma1 * ex1 * phi - gamma1 * sigma * theta1_hat;
    dtheta2 = gamma2 * ex2 * phi - gamma2 * sigma * theta2_hat;

    % Προβολή για a_11 ∈ [-3, -1]
    if theta1_hat(1) <= -3 && dtheta1(1) < 0
        dtheta1(1) = 0;
    elseif theta1_hat(1) >= -1 && dtheta1(1) > 0
        dtheta1(1) = 0;
    end

    % Προβολή για b_2 ∈ [1, ∞)
    if theta2_hat(3) <= 1 && dtheta2(3) < 0
        dtheta2(3) = 0;
    end

   
    theta1_hat = theta1_hat + dt * dtheta1;
    theta2_hat = theta2_hat + dt * dtheta2;

    theta1_hist(:,k+1) = theta1_hat;
    theta2_hist(:,k+1) = theta2_hat;
end

% Plots
figure;
plot(t, x(1,:), 'b', t, x_hat(1,:), '--r'); title('x_1(t) & \hat{x}_1(t)');
legend('x_1', 'x̂_1');
grid on;

figure;
plot(t, x(2,:), 'b', t, x_hat(2,:), '--r'); title('x_2(t) & \hat{x}_2(t)');
legend('x_2', 'x̂_2');
grid on;

figure;
plot(t, x(1,:) - x_hat(1,:), 'g', t, x(2,:) - x_hat(2,:), 'm');
title('e_{x}(t) = x(t) - \hat{x}(t)');
legend('e_{x1}', 'e_{x2}');
grid on;

figure;
plot(t, theta1_hist'); title('\thetâ_1(t)');
legend('â_{11}', 'â_{12}', 'b̂_1');
grid on;

figure;
plot(t, theta2_hist'); title('\thetâ_2(t)');
legend('â_{21}', 'â_{22}', 'b̂_2');
grid on;
