clf; close all

%% Constants
e_0 = 1.60*10^(-19);
epsilon_0 = 8.85*10^(-21);
k_B = 1.38*10^(-23);
d = 100;
T = 300;
epsilon = 80;
c_0 = 0.006;

y_0 = -0.5;
y_d = 0.5;

K = sqrt((2*c_0*(e_0)^2)/(epsilon*epsilon_0*k_B*T));

N = 100;
h = d/(N-1);

%% Matrix A and vector b

% A = zeros(N, N);

% Check the function LU-Decomposition
% A = [10 , -6, 10; 1, -3, 1; 5, -7, 5];
% [L, U] = LU_decomposition(A);
% disp('Matrix L:');
% disp(L);
% disp('Matrix U:');
% disp(U);
% disp('Check:');
% disp(L*U)

y = y_calculation(N, h, y_0, y_d, K);

% Plot c)
x = linspace(0, d, N+2);
analytical_solution = (y_d * sinh(K*x) + y_0 * sinh(K*(d - x))) / sinh(K*d);

figure;
plot(x, y, 'b-', 'LineWidth', 2)
hold on
plot(x, analytical_solution, 'r--', 'LineWidth', 2)
xlabel('x')
ylabel('y(x)')
title('Comparison of analytical and numerical solutions')
legend('Numerical Solution', 'Analytical Solution', Location='best')
grid on

%% (d):
d_val = linspace(1, d, N);

C = zeros(length(d_val), 3);
C_analytical = zeros(length(d_val),3);

dq = zeros(size(d_val));
dy = 1e-7;

for i = 1:N %C(d)
    h_d = d_val(i)/(N-1);

    y_num_d = y_calculation(N, h_d, y_0, y_d, K);
    y_perturbed = y_calculation(N, h_d, y_0 + dy, y_d, K);

    dy_dx_at_0_num = (y_num_d(2) - y_num_d(1)) / h_d;
    dy_dx_at_d_num = (y_num_d(end) - y_num_d(end - 1)) / h_d;

    dy_dx_at_0_perturbed = (y_perturbed(2) - y_perturbed(1)) / h_d;
    dy_dx_at_d_perturbed = (y_perturbed(end) - y_perturbed(end - 1)) / h_d;

    q_num = (2*c_0/(K^2))*(dy_dx_at_d_num - dy_dx_at_0_num);
    q_perturbed = (2*c_0/(K^2))*(dy_dx_at_d_perturbed - dy_dx_at_0_perturbed);
    dq(i) = q_perturbed - q_num;
    
    C(i, 1) = dq(i)/(dy);

    %analytical solution:
    % d0 = (K/sinh(K*d_val(i)))*(y_d - (y_0 + dy)*cosh(K*d_val(i)));
    % dd = (K/sinh(K*d_val(i)))*(y_d*cosh(K*d_val(i)) - (y_0 + dy));
    % q_analytical(i) = (2*c_0)/(K^2) * (dd - d0);
    C_analytical(i, 1) = ((2*c_0)/(K*sinh(K*d_val(i))))*(cosh(d_val(i)*K)-1);
end

k_val = linspace(0.1*K, K, N);
dq_k = zeros(size(d_val));
C_k = zeros(size(d_val));
for k = 1:length(k_val) %C(k)
    y_num_k = y_calculation(N, h, y_0, y_d, k_val(k));
    y_perturbed = y_calculation(N, h, y_0 + dy, y_d, k_val(k));

    dy_dx_at_0_num = (y_num_k(2) - y_num_k(1)) / h;
    dy_dx_at_d_num = (y_num_k(end) - y_num_k(end - 1)) / h;

    dy_dx_at_0_perturbed = (y_perturbed(2) - y_perturbed(1)) / h;
    dy_dx_at_d_perturbed = (y_perturbed(end) - y_perturbed(end - 1)) / h;

    q_num = (2*c_0/((k_val(k))^2))*(dy_dx_at_d_num - dy_dx_at_0_num);
    q_perturbed = (2*c_0/(k_val(k)^2))*(dy_dx_at_d_perturbed - dy_dx_at_0_perturbed);
    dq_k(k) = q_perturbed - q_num;

    C(k, 2) = dq_k(k)/(dy);
    C_analytical(k, 2) = ((2*c_0)/(k_val(k)*sinh(k_val(k)*d)))*(cosh(d*k_val(k))-1);
end

T_val = linspace(1, T, N);
dq_T = zeros(size(d_val));
for t = 1:length(T_val) %C(T)
    k_T = sqrt((2*c_0*(e_0)^2)/(epsilon*epsilon_0*k_B*T_val(t)));

    y_num_k = y_calculation(N, h, y_0, y_d, k_T);
    y_perturbed = y_calculation(N, h, y_0 + dy, y_d, k_T);

    dy_dx_at_0_num = (y_num_k(2) - y_num_k(1)) / h;
    dy_dx_at_d_num = (y_num_k(end) - y_num_k(end - 1)) / h;

    dy_dx_at_0_perturbed = (y_perturbed(2) - y_perturbed(1)) / h;
    dy_dx_at_d_perturbed = (y_perturbed(end) - y_perturbed(end - 1)) / h;

    q_num = (epsilon*epsilon_0*k_B*T_val(t)/(e_0^2))*(dy_dx_at_d_num - dy_dx_at_0_num);
    q_perturbed = (epsilon*epsilon_0*k_B*T_val(t)/(e_0^2))*(dy_dx_at_d_perturbed - dy_dx_at_0_perturbed);
    dq_T(t) = q_perturbed - q_num;

    C(t, 3) = dq_T(t)/(dy);
    C_analytical(t, 3) = ((2*c_0)/(k_T*sinh(k_T*d)))*(cosh(d*k_T)-1);
end

%Plots
figure; hold on
plot(d_val, C(:, 1), 'b-', 'LineWidth', 2)
plot(d_val, C_analytical(:, 1), 'r--', 'LineWidth', 2)
xlabel('Device Width (d, nm)')
ylabel('Capacitance per Unit Area (C)')
title('Capacitance as a Function of Device Width')
legend('Numerical Solution', 'Analytical Solution', Location='best')
grid on


screenSize = get(0, 'ScreenSize');
xPos = (screenSize(3) - 1000) / 2;
yPos = (screenSize(4) - 500) / 2;

figure('Position', [xPos, yPos, 1000, 500]);
subplot(1,2,1)
plot(k_val, C(:, 2), 'b-', 'LineWidth', 2)
hold on
plot(k_val, C_analytical(:, 2), 'r--', 'LineWidth', 2)
xlabel('Debye length (k)')
ylabel('Capacitance per Unit Area (C)')
title('Capacitance as a Function of Debye length ')
legend('Numerical Solution', 'Analytical Solution', Location='best')
grid on

subplot(1,2,2)
plot(T_val, C(:, 3), 'b-', 'LineWidth', 2)
hold on
plot(T_val, C_analytical(:, 3), 'r--', 'LineWidth', 2)
xlabel('Temperatur (T)')
ylabel('Capacitance per Unit Area (C)')
title('Capacitance as a Function of Temperatur')
legend('Numerical Solution', 'Analytical Solution', Location='best')
grid on


%e)

differences = diff(C(:, 1));

threshold = 0.00001;
plateau_i = find(abs(differences) < threshold, 1);

plateau_d = d_val(plateau_i);

fprintf('Optimal d: %.2f\n nm', plateau_d);


%% Functions
function [L, U] = LU_decomposition(A)
    [~, n] = size(A);
    L = eye(n);
    U = A;

    for j = 1:n
        for i = (j+1):n
            L(i, j) = U(i, j)/U(j, j);
            U(i,j:n) = U(i, j:n) - L(i, j)*U(j,j:n);
        end
    end
end

function y = y_calculation(N, h, y_0, y_d, K)
    main_diagonal = (-2/(h^2) - K^2) * ones(N, 1);
    sub_diagonal = 1/(h^2)*ones(N-1, 1);

    A = diag(main_diagonal) + diag(sub_diagonal, 1) + diag(sub_diagonal, -1);

    b = zeros(N, 1);
    %Boundary conditions
    b(1) = -y_0/(h^2);
    b(N) = -y_d/(h^2);
    
    [L, U] = LU_decomposition(A);
    % [L, U] = lu(A);
    %[L, U] = LU_dec_pivoting(A);
    %For N > 150 it makes more sense to use built-in functions. At least on my computer it takes too long
    %y = linsolve(A, b);

    % Solution of the system Lx = b
    x = zeros(1, N);
    x(1) = b(1);
    for i = 2:N
        sum = 0;
        for k = 1:(i-1)
            sum = sum + L(i, k) * x(k);
        end
        x(i) = b(i) - sum;
    end

    % Solution of the system Uy = x
    y = zeros(1, N);
    y(N) = x(N) / U(N, N);
    for i = N-1:-1:1
        sum = 0;
        for k = (i+1):N
            sum = sum + U(i, k) * y(k);
        end
        y(i) = (x(i) - sum) / U(i, i);
    end

    y = horzcat(y_0, y, y_d);
   %for linsolve
   %y = horzcat(y_0, y.', y_d);
end

function E = calculate_energy(C, V)
    E = 0.5 * C * V^2;
end

% function [L, U, P] = LU_dec_pivoting(A)
%     [~, n] = size(A);
%     L = eye(n);
%     P = eye(n);
%     U = A;
% 
%     for j = 1:n
%         [~, pivot_row] = max(abs(U(j:n, j)));
%         pivot_row = pivot_row + j - 1;
% 
%         if pivot_row ~= j
%             % Swap rows in U and P
%             temp = U(j, :);
%             U(j, :) = U(pivot_row, :);
%             U(pivot_row, :) = temp;
% 
%             temp = P(j, :);
%             P(j, :) = P(pivot_row, :);
%             P(pivot_row, :) = temp;
% 
%             if j > 1
%                 % Swap corresponding rows in L
%                 temp = L(j, 1:j-1);
%                 L(j, 1:j-1) = L(pivot_row, 1:j-1);
%                 L(pivot_row, 1:j-1) = temp;
%             end
%         end
% 
%         for i = (j+1):n
%             L(i, j) = U(i, j)/U(j, j);
%             U(i,j:n) = U(i, j:n) - L(i, j)*U(j,j:n);
%         end
%     end
% end