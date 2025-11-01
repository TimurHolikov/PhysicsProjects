%% 1D Schrödinger Equation Solver
close all; clear; clc
tic

%% Parameters
Nx = 800;
Nt = 800;
dx = 1;
dt = 1;

hbar = 1;
mass = 1;

x0 = 20;
sigma = 10;
momentum = 2;

V0_values = [1.0, 1.5, 2.0];
barrier1_pos = 200;
barrier2_pos = 300;
barrier_width = 10;

createGIF = true;
gifDelay = 0.1;

%% Setup
x = linspace(0, (Nx-1)*dx, Nx);
t = 0:dt:Nt*dt;

dx2 = dx^2;
coeff1 = mass * dx2 / hbar^2;
coeff2 = 1i * 2 * mass * dx2 / (hbar * dt);

%% Free Propagation
fprintf('Free propagation...\n');
V_free = zeros(1, Nx);
psi_free = solve_schrodinger(Nx, Nt, V_free, x, sigma, momentum, x0, hbar, coeff1, coeff2);

% Plot free propagation
width_function = @(t) sigma * sqrt(1 + (t / (2 * mass * sigma^2))^2);
expected_velocity = momentum / mass;

screenSize = get(0, 'ScreenSize');
xPos = (screenSize(3) - 700) / 2;
yPos = (screenSize(4) - 500) / 2;
fig1 = figure('Position', [xPos, yPos, 700, 500]);

frameCount = 1;
for ts = 1:5:Nt
    width = width_function(t(ts));
    
    subplot(2, 1, 1)
    plot(x, real(psi_free(:, ts)), 'LineWidth', 1.2)
    grid on; box on
    xlabel('x')
    ylabel('\psi_{real}')
    ylim([min(real(psi_free(:))), max(real(psi_free(:)))])
    xlim([0, 400])
    title(['Time Step: ' num2str(ts) ', Expected velocity: ' num2str(expected_velocity) ', Expected width: ' num2str(width)])
    
    subplot(2, 1, 2)
    plot(x, abs(psi_free(:, ts)).^2, 'r', 'LineWidth', 1.2)
    grid on; box on
    xlabel('x')
    ylabel('|\psi|^2')
    ylim([0, 0.05])
    xlim([0, 400])
    
    drawnow limitrate
    
    if createGIF
        frame = getframe(fig1);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if frameCount == 1
            imwrite(imind, cm, '1D_schrodinger.gif', 'gif', 'DelayTime', gifDelay, 'Loopcount', inf);
        else
            imwrite(imind, cm, '1D_schrodinger.gif', 'gif', 'DelayTime', gifDelay, 'WriteMode', 'append');
        end
        frameCount = frameCount + 1;
    end
end

%% Barrier Propagation
fprintf('Barrier propagation...\n');
[~, V_barrier] = create_potential(V0_values(2), x, barrier1_pos, barrier2_pos, barrier_width);

% Absorbing boundary
capSize = round(0.1 * Nx);
s = linspace(0, 1, capSize);
V_barrier(end-capSize+1:end) = V_barrier(end-capSize+1:end) - 1i * 0.03 * (s.^3);

psi_barrier = solve_schrodinger(Nx, Nt, V_barrier, x, sigma, momentum, x0, hbar, coeff1, coeff2);

% Plot barrier propagation
fig2 = figure('Position', [xPos, yPos, 700, 500]);
frameCount = 1;

for ts = 1:5:Nt
    subplot(2, 1, 1)
    plot(x, real(psi_barrier(:, ts)), 'LineWidth', 1.2)
    hold on
    plot(x, real(V_barrier), 'r', 'LineWidth', 1.5)
    hold off
    grid on; box on
    xlabel('x')
    ylabel('\psi_{real}')
    ylim([-0.2, 0.4])
    xlim([0, 400])
    
    subplot(2, 1, 2)
    plot(x, abs(psi_barrier(:, ts)).^2, 'b', 'LineWidth', 1.2)
    hold on
    plot(x, real(V_barrier) / max(real(V_barrier)) * 0.04, 'r', 'LineWidth', 1.5)
    hold off
    grid on; box on
    xlabel('x')
    ylabel('|\psi|^2')
    ylim([0, 0.04])
    xlim([0, 400])
    
    drawnow limitrate
    
    if createGIF
        frame = getframe(fig2);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if frameCount == 1
            imwrite(imind, cm, '1D_schrodinger_V2_V0=1.5.gif', 'gif', 'DelayTime', gifDelay, 'Loopcount', inf);
        else
            imwrite(imind, cm, '1D_schrodinger_V2_V0=1.5.gif', 'gif', 'DelayTime', gifDelay, 'WriteMode', 'append');
        end
        frameCount = frameCount + 1;
    end
end

fprintf('Total time: %.2f seconds\n', toc);

%% Functions
function psi = solve_schrodinger(Nx, Nt, V, x, sigma, momentum, x0, hbar, coeff1, coeff2)
    psi = complex(zeros(Nx, Nt));
    psi(:, 1) = create_initial_wavepacket(x, sigma, momentum, x0, hbar);
    psi(1, :) = 0;
    psi(Nx, :) = 0;
    
    a = complex(zeros(Nx, 1));
    a(1) = 2 * (1 + coeff1 * V(1) - coeff2);
    for k = 2:Nx-1
        a(k) = 2 * (1 + coeff1 * V(k) - coeff2) - 1/a(k-1);
    end
    
    Omega = complex(zeros(Nx, 1));
    b = complex(zeros(Nx, 1));
    
    for n = 2:Nt
        for k = 2:Nx-1
            Omega(k) = -psi(k-1, n-1) + 2 * (coeff2 + 1 + coeff1 * V(k)) * psi(k, n-1) - psi(k+1, n-1);
        end
        
        b(1) = Omega(1);
        for k = 2:Nx-1
            b(k) = b(k-1) / a(k-1) + Omega(k);
        end
        
        for k = Nx-1:-1:2
            psi(k, n) = (1 / a(k)) * (psi(k+1, n) - b(k));
        end
        
        psi(:, n) = sqrt(1 / trapz(psi(:, n) .* conj(psi(:, n)))) * psi(:, n);
    end
end

function psi0 = create_initial_wavepacket(x, sigma, momentum, x0, hbar)
    prefactor = 1 / (sqrt(sigma * sqrt(pi)));
    exp1 = exp(-(x - x0).^2 / (2 * sigma^2));
    exp2 = exp((1i / hbar) * momentum * x);
    
    psi0 = prefactor * (exp1 .* exp2);
    
    Norm = sqrt(1 / trapz(psi0 .* conj(psi0)));
    psi0 = Norm * psi0;
end

function [V1, V2] = create_potential(V0, x, a, b, d)
    V1 = V0 * (heaviside(x - a) - heaviside(x - (a + d)));
    V2 = V1 + V0 * (heaviside(x - b) - heaviside(x - (b + d)));
end