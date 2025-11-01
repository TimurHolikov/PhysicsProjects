% Assignment 2. Task 3
% Names: Timur Holikov, Lea Ebenschweiger


clf; close all

data = readmatrix('xyzm_dna.txt');
R0 = data(:, 1:3); % coordinates
m = data(:, 4); % masses

Rcut = 5;
k = 1;

N = size(R0, 1);
H = zeros(N, N);

%% (a):
% Hessian matrix
% for i = 1:N
%     for j = i+1:N
%         if (norm(R0(i,:) - R0(j,:))) < Rcut
%             H((3*i-2):(3*i), (3*j-2):(3*j)) = -k * eye(3);
%             H((3*j-2):(3*j), (3*i-2):(3*i)) = -k * eye(3);
%             H((3*i-2):(3*i), (3*i-2):(3*i)) = H((3*i-2):(3*i), (3*i-2):(3*i)) + k * eye(3);
%             H((3*j-2):(3*j), (3*j-2):(3*j)) = H((3*j-2):(3*j), (3*j-2):(3*j)) + k * eye(3);
%         end
%     end
% end

for i = 1:N
    for j = 1:N
        if i == j
            H(i, j) = k;
        elseif norm(R0(i, :) - R0(j, :)) < Rcut
            H(i, j) = -k;
        end
    end
end

%% (b):
% M = diag(repmat(m, 3, 1));
M = diag(m);
K = M^(-1/2) * H * M^(-1/2);

eigenvalues = zeros(10,1);
EVs = zeros(size(K, 1), 10);

%Check if the matrix K diagonalizable
    % Compute the eigenvalues and eigenvectors of the matrix
    [V, D] = eig(K);
    
    % Check if the matrix V of eigenvectors has linearly independent columns
    if rank(V) == size(K, 1)
        disp('The matrix is diagonalizable.');
    else
        error('The matrix is not diagonalizable. The power method is not applicable');
    end

K_Def = K;
for i = 1:10
    x = rand(size(K_Def, 1), 1); % Initial vector
    x = x / norm(x); % Normalize

    v = power_method(K_Def,x);

    lambda = (v' * K_Def * v)/(v' * v);
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

% [VV, DD] = eigs(K,10)
%% (c)

% EVs_orth = orth(EVs);
% EVs_norm = normalize(EVs_orth);

% EVs_reshaped = reshape(EVs, N, 3, 10);
% 
% figure
% for i = 1:10
%     subplot(5, 2, i)
%     plot(R0(:, 3), EVs_reshaped(:, 3, i));
%     title(['EV for \lambda = ', num2str(eigenvalues(i))]);
%     xlabel('z');
%     ylabel('Eigenvector');
% end

screenSize = get(0, 'ScreenSize');
xPos = (screenSize(3) - 1300) / 2;
yPos = (screenSize(4) - 800) / 2;

figure('Position', [xPos, yPos, 1300, 800])
for i = 1:10
    subplot(5, 2, i)
    plot(R0(:, 3), EVs(:,i), 'LineWidth', 2)
    title(['EV for \lambda = ', num2str(eigenvalues(i))])
    xlabel('z')
    ylabel('Eigenvector')
end

%% (d)
n = size(K, 1);

row_sums = sum(abs(K), 2);
diag_elements = diag(K);

ri = row_sums - abs(diag_elements);

min_eig = min(diag_elements - ri);
max_eig = max(diag_elements + ri);

disp(['Minimum Eigenvalue: ', num2str(min_eig)]);
disp(['Maximum Eigenvalue: ', num2str(max_eig)]);

% Gershgorin disks for rows
% for i = 1:n
%     ri = sum(abs(K(i, :))) - abs(K(i, i));
%     min_eig = min(min_eig, ri - K(i,i));
%     max_eig = max(min_eig, ri + K(i,i));
% end

%% (e):
K_prime = construct_matrix(K);
eig_K_prime = eig(K_prime);

if all(abs(eig_K_prime) < 1)
    disp('All absolute values of the eigenvalues are strictly smaller than 1');
else
    disp('Not all absolute values of the eigenvalues are strictly smaller than 1');
end

disp(['The smallest eigenvalue of K_prime: ', num2str(min(eig_K_prime))]);
%% (f):
K_p_I = K_prime + eye(size(K_prime));

K_p_I_inv = inv(K_p_I);

% Apply the power method to (K_prime + I)^-1
x_p = rand(size(K_prime, 1), 1);
x_p = x_p / norm(x_p);

u = power_method(K_p_I_inv,x_p);

lambda_K_p_I = (u' * K_p_I_inv * u)/(u' * u);

disp(['Eigenvalue of (K_prime + I)^-1: ', num2str(lambda_K_p_I)]);
disp(['The smallest eigenvalue of K_prime = 1/(EV((K_prime + I)^-1)) - 1 = ', num2str(min(inv(lambda_K_p_I) - 1))]);

%% (e): The Neumann series
I_t = eye(size(K_prime));

for p = 1:1e3
    I_t = I_t + ((-1*K_prime)^p)*(K_p_I);
end

disp('For p -> inf: I_t -> I ');

K_I_inv = eye(size(K_prime)); 
for p = 1:1e4
    K_I_inv = K_I_inv + ((-1*K_prime)^p);
    if max(max(abs((K_prime)^p))) < 1e-8
        break;
    end
end

if (norm(K_I_inv - K_p_I_inv)) < 1e-4
    disp('(K_prime + I)^-1 computed with Neumann series is equal to built-in implementation')
end

%% (h):
K_p_I_inv_Def = K_p_I_inv;

lambda_K_prime = zeros(10,1);
s_eigenvalues_K_prime = zeros(10,1);
eigenvectors_K_prime = zeros(size(K_p_I_inv, 1), 10);

% Smallest value

for i = 1:10
    % Apply the power method to (K_prime + I)^-1
    x_rand = rand(size(K_prime, 1), 1);
    x_rand = x_rand / norm(x_rand);

    uu = power_method(K_p_I_inv_Def,x_rand);

    lambda_K_prime = (uu' * K_p_I_inv_Def * uu)/(uu' * uu);

    s_eigenvalues_K_prime(i) = 1/(lambda_K_prime) - 1;

    % Gram-Schmidt process
    for j = 1:i-1
        uu = uu - (uu' *  eigenvectors_K_prime(:, j)) *  eigenvectors_K_prime(:, j);
    end
    uu = uu / norm(uu);

    eigenvectors_K_prime(:, i) = uu;

    % Deflation
    K_p_I_inv_Def = K_p_I_inv_Def - lambda_K_prime * (uu * uu');
end

% Smallest magnitude:

% K_prime_Def = K_prime;
% sigma = 1e-5; % shift close to zero
% 
% lambda_K_prime = zeros(10,1);
% s_eigenvalues_K_prime = zeros(10,1);
% eigenvectors_K_prime = zeros(size(K_prime, 1), 10);
% 
% %LU Decomposition (K_prime - sigma*I)
% [L, U] = lu(K_prime - sigma*eye(size(K_prime)));
% 
% for i = 1:10
%     x_rand = rand(size(K_prime, 1), 1);
%     x_rand = x_rand / norm(x_rand);
% 
%     for j = 1:1000 % number of iterations
%         y = U \ (L \ x_rand); % (K_prime - sigma*I)y = x_rand
%         x_rand = y / norm(y); % normalize
%     end
% 
%     uu = x_rand;
%     lambda_K_prime = (uu' * K_prime * uu)/(uu' * uu);
% 
%     s_eigenvalues_K_prime(i) = lambda_K_prime;
% 
%     % Gram-Schmidt process
%     for j = 1:i-1
%         uu = uu - (uu' *  eigenvectors_K_prime(:, j)) *  eigenvectors_K_prime(:, j);
%     end
%     uu = uu / norm(uu);
% 
%     eigenvectors_K_prime(:, i) = uu;
% 
%     % Deflation
%     K_prime_Def = K_prime_Def - lambda_K_prime * (uu * uu');
% end

%% (i):
figure('Position', [xPos, yPos, 1300, 800])
for i = 1:10
    if eigenvectors_K_prime(i) ~= 0
        subplot(5,2,i)
        plot(R0(1:N/2, 3), eigenvectors_K_prime(1:N/2, i), 'Color', 'r') % first DNA strand
        hold on

        plot(R0(N/2+1:end, 3), eigenvectors_K_prime(N/2+1:end, i), 'Color', 'b') % second DNA strand
        title(['Eigenvector for \lambda = ', num2str(s_eigenvalues_K_prime(i))]);
        xlabel('z');
        ylabel('Eigenvector');
        hold off
    end
end

%% Functions:
function x_final = power_method(A,x)
    x_final = x;

    for n = 1:1e5
        x_old=x_final;

        x_final=A*x_old;
        x_final = x_final/norm(x_final);

        if norm(x_final-x_old) < eps()
            break;
        end
    end
end

function K_prime = construct_matrix(K)
    max_singular_value = norm(K, 2);
    scaling_factor = max_singular_value * 1.01; % Define a scaling factor slightly larger than the maximum absolute eigenvalue

    K_prime = K / scaling_factor;
end
