% Assignment 2. Task 4
% Names: Timur Holikov, Lea Ebenschweiger


% Assignment 2. Task 4
% Names: Timur Holikov, Lea Ebenschweiger


clf; close all

data = readmatrix('xyzm_dna.txt');
R0 = data(:, 1:3); % coordinates
m = data(:, 4); % masses

Rcut = 5;
k = 1;

N = size(R0, 1);
H = zeros(3*N, 3*N);

%% (a):
% Hessian matrix
for i = 1:N
    for j = i+1:N
        if (norm(R0(i,:) - R0(j,:))) < Rcut
            H((3*i-2):(3*i), (3*j-2):(3*j)) = -k * eye(3);
            H((3*j-2):(3*j), (3*i-2):(3*i)) = -k * eye(3);
            H((3*i-2):(3*i), (3*i-2):(3*i)) = H((3*i-2):(3*i), (3*i-2):(3*i)) + k * eye(3);
            H((3*j-2):(3*j), (3*j-2):(3*j)) = H((3*j-2):(3*j), (3*j-2):(3*j)) + k * eye(3);
        end
    end
end

% Eigenvalues of K:
M = diag(repmat(m, 3, 1));
K = M^(-1/2) * H * M^(-1/2);

eigenvalues = zeros(3,1);
eigenvectors = zeros(size(K, 1), 3);

% 3 largest
for i = 1:3
    x = rand(size(K, 1), 1); % Initial vector
    x = x / norm(x); % Normalize

    v = power_method(K,x);

    lambda = (v' * K * v)/(v' * v);
    eigenvalues(i) = lambda;

    % Gram-Schmidt process
    for j = 1:i-1
        v = v - (v' * eigenvectors(:, j)) * eigenvectors(:, j);
    end
    v = v / norm(v);

    eigenvectors(:, i) = v;

    % Deflation
    K = K - lambda * (v * v');
end

%3 smallest
K_prime = construct_matrix(K);
K_inv = inv(K_prime+eye(size(K)));

eigenvalues_small = zeros(3,1);
eigenvectors_small = zeros(size(K, 1), 3);

for i = 1:3
    xx = rand(size(K_inv, 1), 1); % Initial vector
    xx = xx / norm(xx); % Normalize

    u = power_method(K_inv,xx);

    lambda_prime = (u' * K_inv * u)/(u' * u);

    % Gram-Schmidt process
    for j = 1:i-1
        u = u - (u' * eigenvectors_small(:, j)) * eigenvectors_small(:, j);
    end
    u = u / norm(u);

    eigenvalues_small(i) = 1/lambda_prime - 1;

    eigenvectors_small(:, i) = u;

    % Deflation
    K_inv = K_inv - lambda_prime * (u * u');
end


%% (b):

% % Superposition of all modes
% t = linspace(0, 2*pi, 400);
% 
% frames = cell(length(t), 1);
% 
% filename = 'animation_superposition.gif';
% h = figure;
% axis tight manual
% set(gca,'nextplot','replacechildren','visible','off')
% 
% % Superposition of all the eigenmodes. Commented, because takes too much time to calculate. Results were exported as gif files
% % while true
%     for i = 1:length(t)
%         displacement = zeros(size(R0));
% 
%         % Loop over all modes
%         for mode = 1:size(eigenvectors, 2)
%             % Calculate displacement for current mode
%             displacement_mode = sin(t(i) + mode) * reshape(eigenvectors(:, mode), 3, [])';
% 
%             displacement = displacement + displacement_mode;
%         end
% 
%         R = R0 + displacement;
%         clf;
% 
%         plot3(R(1:N/2, 1), R(1:N/2, 2), R(1:N/2, 3), 'ro'); % First DNA strand
%         hold on
%         plot3(R(N/2+1:end, 1), R(N/2+1:end, 2), R(N/2+1:end, 3), 'bo'); % Second DNA
%         plot3(R0(:, 1), R0(:, 2), R0(:, 3), 'k'); % Equilibrium position
%         % grid on
% 
%         xlim([min(R0(:, 1)), max(R0(:, 1))]);
%         ylim([min(R0(:, 2)), max(R0(:, 2))]);
%         zlim([min(R0(:, 3)), max(R0(:, 3))]);
% 
%             view([-48.0 11.0])
% 
%         % Capture the plot as an image
%         frame = getframe(h);
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
% 
%         % Store the frame in the cell array
%         frames{i} = imind;
%     end
% % end
% 
% % Write all the frames to the GIF file at once
% for idx = 1:length(frames)
%     if idx == 1
%         imwrite(frames{idx},cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.00001);
%     else
%         imwrite(frames{idx},cm,filename,'gif','WriteMode','append', 'DelayTime', 0.00001);
%     end
% end
% 
% close(h);

% First mode ?
% for i = 1:length(t)
%     % Calculate new positions
%     R = R0 + sin(t(i)) * reshape(eigenvectors(:, 1), 3, [])';
% 
%     clf;
% 
%     plot3(R(:, 1), R(:, 2), R(:, 3), 'o');
% 
%     % Set axis limits
%     xlim([min(R0(:, 1)), max(R0(:, 1))]);
%     ylim([min(R0(:, 2)), max(R0(:, 2))]);
%     zlim([min(R0(:, 3)), max(R0(:, 3))]);
% 
%     % Pause for a short time to create animation
%     pause(0.01);
% end

t = linspace(0, 2*pi, 1000); % Increase the number of points for smoother animation

% Create figures
h_large = figure;
h_small = figure;

% Loop over all time points
for i = 1:length(t)
    % Loop over all modes
    for mode = 1:3
        % Calculate displacement for current mode
        displacement_mode_large = sin(t(i) + mode) * reshape(eigenvectors(:, mode), 3, [])';
        displacement_mode_small = sin(t(i) + mode) * reshape(eigenvectors_small(:, mode), 3, [])';

        % Calculate new positions
        R_large = R0 + displacement_mode_large;
        R_small = R0 + displacement_mode_small;

        % Plot large modes
        figure(h_large);
        subplot(1,3,mode);
        plot3(R_large(1:N/2, 1), R_large(1:N/2, 2), R_large(1:N/2, 3), 'ro'); % First DNA strand
        hold on
        plot3(R_large(N/2+1:end, 1), R_large(N/2+1:end, 2), R_large(N/2+1:end, 3), 'bo'); % Second DNA strand
        plot3(R0(:, 1), R0(:, 2), R0(:, 3), 'k'); % Equilibrium position
        title(['Large Eigenmode ', num2str(mode)]);

        % Plot small modes
        figure(h_small);
        subplot(1,3,mode);
        plot3(R_small(1:N/2, 1), R_small(1:N/2, 2), R_small(1:N/2, 3), 'ro'); % First DNA strand
        hold on
        plot3(R_small(N/2+1:end, 1), R_small(N/2+1:end, 2), R_small(N/2+1:end, 3), 'bo'); % Second DNA strand
        plot3(R0(:, 1), R0(:, 2), R0(:, 3), 'k'); % Equilibrium position
        title(['Small Eigenmode ', num2str(mode)]);
        
        drawnow;

        pause(0.01); % Pause for a short time to create animation
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
