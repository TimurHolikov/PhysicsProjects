% DNA Structure Eigenvalue Analysis
clf; close all;

%% Load Data
data = readmatrix('xyzm_dna.txt');
R0 = data(:, 1:3); % coordinates
m = data(:, 4); % masses

Rcut = 5;
k = 1;
N = size(R0, 1);

%% (a): Hessian matrix
fprintf('Building Hessian matrix...\n');
tic;

% Vectorized distance calculation
dist_matrix = pdist2(R0, R0);
H = zeros(N, N);

for i = 1:N
    for j = 1:N
        if i == j
            H(i, j) = k;
        elseif dist_matrix(i, j) < Rcut
            H(i, j) = -k;
        end
    end
end

fprintf('Hessian built in %.2f seconds\n', toc);

%% (b): Eigenvalue computation
M = diag(m);
K = M^(-1/2) * H * M^(-1/2);

eigenvalues = zeros(10, 1);
EVs = zeros(size(K, 1), 10);

% Check if the matrix K diagonalizable
[V, D] = eig(K);

if rank(V) == size(K, 1)
    disp('The matrix is diagonalizable.');
else
    error('The matrix is not diagonalizable. The power method is not applicable');
end

K_Def = K;
for i = 1:10
    x = rand(size(K_Def, 1), 1); % Initial vector
    x = x / norm(x); % Normalize
    
    v = power_method(K_Def, x);
    
    lambda = (v' * K_Def * v) / (v' * v);
    eigenvalues(i) = lambda;
    
    % Gram-Schmidt process
    for j = 1:i-1
        v = v - (v' * EVs(:, j)) * EVs(:, j);
    end
    v = v / norm(v);
    
    EVs(:, i) = v;
    
    % Deflation
    K_Def = K_Def - lambda * (v * v');
end

%% (c): Plot eigenvectors
screenSize = get(0, 'ScreenSize');
xPos = (screenSize(3) - 1300) / 2;
yPos = (screenSize(4) - 800) / 2;

figure('Position', [xPos, yPos, 1300, 800])
for i = 1:10
    subplot(5, 2, i)
    plot(R0(:, 3), EVs(:, i), 'LineWidth', 2)
    title(['EV for \lambda = ', num2str(eigenvalues(i))])
    xlabel('z')
    ylabel('Eigenvector')
end

%% (d): Gershgorin bounds
n = size(K, 1);

row_sums = sum(abs(K), 2);
diag_elements = diag(K);

ri = row_sums - abs(diag_elements);

min_eig = min(diag_elements - ri);
max_eig = max(diag_elements + ri);

disp(['Minimum Eigenvalue: ', num2str(min_eig)]);
disp(['Maximum Eigenvalue: ', num2str(max_eig)]);

%% (e): Construct K_prime
K_prime = construct_matrix(K);
eig_K_prime = eig(K_prime);

if all(abs(eig_K_prime) < 1)
    disp('All absolute values of the eigenvalues are strictly smaller than 1');
else
    disp('Not all absolute values of the eigenvalues are strictly smaller than 1');
end

disp(['The smallest eigenvalue of K_prime: ', num2str(min(eig_K_prime))]);

%% (f): Inverse power method
K_p_I = K_prime + eye(size(K_prime));

K_p_I_inv = inv(K_p_I);

% Apply the power method to (K_prime + I)^-1
x_p = rand(size(K_prime, 1), 1);
x_p = x_p / norm(x_p);

u = power_method(K_p_I_inv, x_p);

lambda_K_p_I = (u' * K_p_I_inv * u) / (u' * u);

disp(['Eigenvalue of (K_prime + I)^-1: ', num2str(lambda_K_p_I)]);
disp(['The smallest eigenvalue of K_prime = 1/(EV((K_prime + I)^-1)) - 1 = ', num2str(min(inv(lambda_K_p_I) - 1))]);

%% (g): The Neumann series
fprintf('Computing Neumann series...\n');
tic;

K_I_inv = eye(size(K_prime));
K_power = eye(size(K_prime));
K_neg = -K_prime;

for p = 1:1e4
    K_power = K_power * K_neg;
    K_I_inv = K_I_inv + K_power;
    if max(max(abs(K_power))) < 1e-8
        break;
    end
end

if (norm(K_I_inv - K_p_I_inv)) < 1e-4
    disp('(K_prime + I)^-1 computed with Neumann series is equal to built-in implementation')
end

fprintf('Neumann series completed in %.2f seconds\n', toc);

%% (h): Compute smallest eigenvalues
K_p_I_inv_Def = K_p_I_inv;

lambda_K_prime = zeros(10, 1);
s_eigenvalues_K_prime = zeros(10, 1);
eigenvectors_K_prime = zeros(size(K_p_I_inv, 1), 10);

for i = 1:10
    % Apply the power method to (K_prime + I)^-1
    x_rand = rand(size(K_prime, 1), 1);
    x_rand = x_rand / norm(x_rand);
    
    uu = power_method(K_p_I_inv_Def, x_rand);
    
    lambda_K_prime(i) = (uu' * K_p_I_inv_Def * uu) / (uu' * uu);
    
    s_eigenvalues_K_prime(i) = 1 / (lambda_K_prime(i)) - 1;
    
    % Gram-Schmidt process
    for j = 1:i-1
        uu = uu - (uu' * eigenvectors_K_prime(:, j)) * eigenvectors_K_prime(:, j);
    end
    uu = uu / norm(uu);
    
    eigenvectors_K_prime(:, i) = uu;
    
    % Deflation
    K_p_I_inv_Def = K_p_I_inv_Def - lambda_K_prime(i) * (uu * uu');
end

%% (i): Plot DNA strand eigenvectors
figure('Position', [xPos, yPos, 1300, 800])
for i = 1:10
    if eigenvectors_K_prime(i) ~= 0
        subplot(5, 2, i)
        plot(R0(1:N/2, 3), eigenvectors_K_prime(1:N/2, i), 'Color', 'r') % first DNA strand
        hold on
        
        plot(R0(N/2+1:end, 3), eigenvectors_K_prime(N/2+1:end, i), 'Color', 'b') % second DNA strand
        title(['Eigenvector for \lambda = ', num2str(s_eigenvalues_K_prime(i))]);
        xlabel('z');
        ylabel('Eigenvector');
        hold off
    end
end

fprintf('\nAnalysis complete!\n');

%% Functions
function x_final = power_method(A, x)
    x_final = x;
    
    for n = 1:1e5
        x_old = x_final;
        
        x_final = A * x_old;
        x_final = x_final / norm(x_final);
        
        if norm(x_final - x_old) < eps()
            break;
        end
    end
end

function K_prime = construct_matrix(K)
    max_singular_value = norm(K, 2);
    scaling_factor = max_singular_value * 1.01;
    K_prime = K / scaling_factor;
end