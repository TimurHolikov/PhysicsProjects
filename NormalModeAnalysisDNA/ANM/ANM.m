%% DNA Normal Mode Analysis (Anisotropic)
clf; close all; clear; clc;
fprintf('--- DNA Anisotropic Normal Mode Analysis ---\n');

%% Load Data
data = readmatrix('xyzm_dna.txt');
R0 = data(:, 1:3);
m  = data(:, 4);
N  = size(R0, 1);

Rcut = 5; k = 1;

%% Build anisotropic Hessian
fprintf('Building anisotropic Hessian...\n');
tic;
H = zeros(3*N);

for i = 1:N
    for j = i+1:N
        rij = R0(j,:) - R0(i,:);
        dist = norm(rij);
        if dist < Rcut && dist > 0
            u = rij / dist;
            Hblock = -k * (u' * u);
            idx_i = (3*i-2):(3*i);
            idx_j = (3*j-2):(3*j);
            H(idx_i, idx_j) = H(idx_i, idx_j) + Hblock;
            H(idx_j, idx_i) = H(idx_j, idx_i) + Hblock;
            H(idx_i, idx_i) = H(idx_i, idx_i) - Hblock;
            H(idx_j, idx_j) = H(idx_j, idx_j) - Hblock;
        end
    end
end
fprintf('Hessian built in %.2f s\n', toc);

%% Build mass matrix
M = diag(repmat(m, 3, 1));

%% Solve generalized eigenproblem H*phi = omega^2 * M * phi
fprintf('Solving generalized eigenproblem...\n');
tic;
opts.issym = true; opts.isreal = true;

[phi_small, Om2_small] = eigs(H, M, 3, 'smallestreal', opts);
[phi_large, Om2_large] = eigs(H, M, 3, 'largestreal', opts);

omega_small = sqrt(max(diag(Om2_small), 0));
omega_large = sqrt(max(diag(Om2_large), 0));

fprintf('Eigenproblem solved in %.2f s\n', toc);
fprintf('Smallest ω: %.4f, %.4f, %.4f\n', omega_small);
fprintf('Largest  ω: %.4f, %.4f, %.4f\n', omega_large);

%% Normalize for visualization
phi_small = phi_small ./ max(abs(phi_small(:)));
phi_large = phi_large ./ max(abs(phi_large(:)));

%% Animation parameters
t = linspace(0, 2*pi, 100)
scale_small = 15;
scale_large = 8;

xlim_ = [min(R0(:,1))-1, max(R0(:,1))+1];
ylim_ = [min(R0(:,2))-1, max(R0(:,2))+1];
zlim_ = [min(R0(:,3))-1, max(R0(:,3))+1];

%% Functions
function plot_dna(R, R0, N)
    plot3(R0(:,1), R0(:,2), R0(:,3), '-', 'Color',[0.8 0.8 0.8], 'LineWidth',0.6); % equilibrium
    plot3(R(1:N/2,1), R(1:N/2,2), R(1:N/2,3), '-', 'Color',[1 0.25 0.25], 'LineWidth',1.4);
    plot3(R(N/2+1:end,1), R(N/2+1:end,2), R(N/2+1:end,3), '-', 'Color',[0.25 0.55 1], 'LineWidth',1.4);
    scatter3(R(1:N/2,1), R(1:N/2,2), R(1:N/2,3), 10, [1 0.3 0.3], 'filled');
    scatter3(R(N/2+1:end,1), R(N/2+1:end,2), R(N/2+1:end,3), 10, [0.3 0.5 1], 'filled');
end

% Animate smallest eigenmodes
fprintf('Animating smallest modes...\n');
figure('Color',[0.08 0.08 0.08],'Position',[100 100 1400 400]);
gif_small = 'small_eigenmodes_animation.gif';

for f = 1:length(t)
    for mode = 1:3
        phi = reshape(phi_small(:, mode), 3, []).';
        R = R0 + scale_small * sin(t(f)) * phi;

        subplot(1,3,mode); cla; hold on;
        plot_dna(R, R0, N);
        title(sprintf('Small Mode %d (ω = %.4f)', mode, omega_small(mode)), 'Color','w','FontWeight','bold');
        xlabel('X','Color','w'); ylabel('Y','Color','w'); zlabel('Z','Color','w');
        xlim(xlim_); ylim(ylim_); zlim(zlim_);
        grid on; axis vis3d;
        set(gca,'Color',[0.08 0.08 0.08],'XColor','w','YColor','w','ZColor','w');
        view([-48 11]);
    end
    drawnow;

    frame = getframe(gcf);
    [im, map] = rgb2ind(frame2im(frame), 256);
    if f == 1
        imwrite(im, map, gif_small, 'gif', 'Loopcount', inf, 'DelayTime', 0.04);
    else
        imwrite(im, map, gif_small, 'gif', 'WriteMode', 'append', 'DelayTime', 0.04);
    end
end
fprintf('Saved %s\n', gif_small);

% Animate largest eigenmodes
fprintf('Animating largest modes...\n');
figure('Color',[0.08 0.08 0.08],'Position',[100 600 1400 400]);
gif_large = 'large_eigenmodes_animation.gif';

for f = 1:length(t)
    for mode = 1:3
        phi = reshape(phi_large(:, mode), 3, []).';
        R = R0 + scale_large * sin(t(f)) * phi;

        subplot(1,3,mode); cla; hold on;
        plot_dna(R, R0, N);
        title(sprintf('Large Mode %d (ω = %.4f)', mode, omega_large(mode)), 'Color','w','FontWeight','bold');
        xlabel('X','Color','w'); ylabel('Y','Color','w'); zlabel('Z','Color','w');
        xlim(xlim_); ylim(ylim_); zlim(zlim_);
        grid on; axis vis3d;
        set(gca,'Color',[0.08 0.08 0.08],'XColor','w','YColor','w','ZColor','w');
        view([-48 11]);
    end
    drawnow;

    frame = getframe(gcf);
    [im, map] = rgb2ind(frame2im(frame), 256);
    if f == 1
        imwrite(im, map, gif_large, 'gif', 'Loopcount', inf, 'DelayTime', 0.04);
    else
        imwrite(im, map, gif_large, 'gif', 'WriteMode', 'append', 'DelayTime', 0.04);
    end
end
fprintf('Saved %s\n', gif_large);
