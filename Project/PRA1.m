clear; clc;

% Χρονικά διαστήματα
T = 20;              % διάρκεια προσομοίωσης
dt = 0.01;           % χρονικό βήμα
t = 0:dt:T;
N = length(t);

% Είσοδος u(t)
u = sin(2*t);      % εναλλασσόμενο σήμα εισόδου

% Πραγματικά A, B όπως δίνει η εκφώνηση
A_true = [-2.15  0.25;
          -0.75 -2.00];

B_true = [0;
          1.5];

% Θέσεις αληθινών παραμέτρων
theta1_star = [A_true(1,1); A_true(1,2); B_true(1)];
theta2_star = [A_true(2,1); A_true(2,2); B_true(2)];

% Αρχικές καταστάσεις
x = zeros(2, N);           % πραγματικό x
x(:,1) = [1; 1];           % αρχική συνθήκη

x_hat = zeros(2, N);       % εκτιμημένο x
x_hat(:,1) = [0; 0];       % αρχική τιμή εκτιμητή

% Εκτιμητές παραμέτρων
theta1_hat = [-2; 0; 0];   % a11 έχει περιορισμό [-3, -1]
theta2_hat = [0; 0; 1];    % b2 έχει περιορισμό [1, ∞)

% Αποθήκευση ιστορίας εκτιμήσεων
theta1_hist = zeros(3,N);
theta2_hist = zeros(3,N);
theta1_hist(:,1) = theta1_hat;
theta2_hist(:,1) = theta2_hat;

% Learning rates και stabilizing terms
gamma1 = 5;
gamma2 = 5;
theta_m1 = 3;
theta_m2 = 3;

% Πίνακες για σφάλματα
e_x = zeros(2,N);

% Προσομοίωση
for k = 1:N-1
    phi = [x(1,k); x(2,k); u(k)];

    % Σύστημα
    x_dot = A_true * x(:,k) + B_true * u(k);
    x(:,k+1) = x(:,k) + dt * x_dot;

    % Εκτιμητής κατάστασης
    e_x_k = x(:,k) - x_hat(:,k);
    x_hat_dot = [phi' * theta1_hat + theta_m1 * e_x_k(1);
                 phi' * theta2_hat + theta_m2 * e_x_k(2)];
    x_hat(:,k+1) = x_hat(:,k) + dt * x_hat_dot;

    % Προσαρμοστικός νόμος για θ1 (a11 έχει περιορισμό)
    dtheta1 = gamma1 * e_x_k(1) * phi;

    if theta1_hat(1) > -3 && theta1_hat(1) < -1
        % μέσα στο εσωτερικό του διαστήματος
        dtheta1(1) = dtheta1(1);
    elseif theta1_hat(1) <= -3 && dtheta1(1) > 0
        dtheta1(1) = dtheta1(1);
    elseif theta1_hat(1) >= -1 && dtheta1(1) < 0
        dtheta1(1) = dtheta1(1);
    else
        dtheta1(1) = 0;
    end

    theta1_hat = theta1_hat + dt * dtheta1;

    % Προσαρμοστικός νόμος για θ2 (b2 έχει περιορισμό)
    dtheta2 = gamma2 * e_x_k(2) * phi;

    if theta2_hat(3) > 1
        dtheta2(3) = dtheta2(3);
    elseif theta2_hat(3) == 1 && dtheta2(3) >= 0
        dtheta2(3) = dtheta2(3);
    else
        dtheta2(3) = 0;
    end

    theta2_hat = theta2_hat + dt * dtheta2;

    % Αποθήκευση
    theta1_hist(:,k+1) = theta1_hat;
    theta2_hist(:,k+1) = theta2_hat;
    e_x(:,k) = e_x_k;
end


figure;
plot(t, x(1,:), 'b', t, x_hat(1,:), 'r--');
legend('x_1', 'x̂_1'); title('Κατάσταση x_1 και εκτίμηση'); grid on;

figure;
plot(t, x(2,:), 'b', t, x_hat(2,:), 'r--');
legend('x_2', 'x̂_2'); title('Κατάσταση x_2 και εκτίμηση'); grid on;

figure;
plot(t, theta1_hist(1,:), t, theta1_hist(2,:), t, theta1_hist(3,:));
legend('â_{11}', 'â_{12}', 'b̂_1'); title('Εκτιμήσεις θ̂₁');

figure;
plot(t, theta2_hist(1,:), t, theta2_hist(2,:), t, theta2_hist(3,:));
legend('â_{21}', 'â_{22}', 'b̂_2'); title('Εκτιμήσεις θ̂₂');

figure;
plot(t, e_x(1,:), t, e_x(2,:));
legend('e_{x1}', 'e_{x2}'); title('Σφάλματα κατάστασης'); grid on;

