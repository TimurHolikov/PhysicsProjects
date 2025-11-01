close all; clear; clc
tic
%% Options
flagCreateGIF = false;   % set to true to write GIFs (slower)
gifDelay = 0.1;          % delay between GIF frames

flagA1 = 0; % set to 1 to enable GIF creation
delay1 = 0.1;
cs = 1;

%% Parameters
Nx = 800;        % Number of spatial grid points
Nt = 800;       % Number of time steps       
h = 1;
m = 1;

x_0 = 0;
sigma = 10;
q = 2;

V_0 = [1.0, 1.5, 2.0];
aa = 200;
bb = 300;
dd = 10;

dx = 1; x = linspace(-50, (Nx-1)*dx, Nx);
dt = 1; t = 0:dt:Nt*dt;

%% Free Propagation
V = complex(zeros(1, Nx));

[psi, Omega] = crank_nicolson_solve(Nx, dx, Nt, dt, V, m, x, sigma, q, x_0, h);

%% Velocity and width
width = zeros(1, Nt);
width_function = @(sigma, tt, m) sigma * sqrt(1 + (tt / (2 * m * sigma^2))^2);

%% Plots

screenSize = get(0, 'ScreenSize');
xPos = (screenSize(3) - 700) / 2;
yPos = (screenSize(4) - 500) / 2;

figure('Position', [xPos, yPos, 700, 500]);
for ts = 1:9:Nt
    width(ts) = width_function(sigma, ts, m);

    subplot(2, 1, 1)
    plot(x, real(psi(:, ts)))
    grid on
    box on
    xlabel('x')
    ylabel('\psi_{real}')
    ylim([min(real(psi(:))), max(real(psi(:)))])
    xlim([0,400])
    title(['Time Step: ' num2str(ts) ', Expected velocity: ' num2str(q/m) ', Expected width: ' num2str(width(ts))])

    subplot(2, 1, 2)
    plot(x, abs(psi(:, ts)).^2,'r')
    grid on
    box on
    xlabel('x')
    ylabel('|\psi|^2')
    ylim([0, 0.05])
    xlim([0,400])

    pause(0.1)

    if flagA1 == 1
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % On the first loop, create the file. In subsequent loops, append
        if cs == 1
            imwrite(imind,cm,'1D_schrodinger.gif','gif','DelayTime',delay1,'loopcount',inf);
        else
            imwrite(imind,cm,'1D_schrodinger.gif','gif','DelayTime',delay1,'writemode','append');
        end
        cs = cs+1;
    end
end


%% Barriers
[V1, V2] = potential(V_0(2), x, aa, bb, dd);

%% Crank-Nicolson solver
[psi_b, ~] = crank_nicolson_solve(Nx, dx, Nt, dt, V2, m, x, sigma, q, x_0, h);

figure
for ts = 1:9:Nt

    subplot(2, 1, 1)
    plot(x, real(psi_b(:, ts)))
    hold on
    plot(x, V2, 'r')
    hold off
    grid on
    box on
    xlabel('x')
    ylabel('\psi_{real}')
    ylim([-0.2, 0.4])
    xlim([0, 400])

    subplot(2, 1, 2)
    plot(x, abs(psi_b(:, ts)).^2,'b')
    hold on
    plot(x, V2, 'r')
    hold off
    grid on
    box on
    xlabel('x')
    ylabel('|\psi|^2')
    ylim([0, 0.04])
    xlim([0,400])

    pause(0.1)

    if flagA1 == 1
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % On the first loop, create the file. In subsequent loops, append
        if cs == 1
            imwrite(imind,cm,'1D_schrodinger_V2_V0=1.5.gif','gif','DelayTime',delay1,'loopcount',inf);
        else
            imwrite(imind,cm,'1D_schrodinger_V2_V0=1.5.gif','gif','DelayTime',delay1,'writemode','append');
        end
        cs = cs+1;
    end
end

toc

%% Functions
function [psi, Omega] = crank_nicolson_solve(Nx, dx, Nt, dt, V, m, x, sigma, q, x_0, h)
    psi = complex(zeros(Nx, Nt));
    psi(:, 1) = initial_psi0(x, sigma, q, x_0, h);
    psi(1, :) = 0; psi(Nx, :) = 0; % Boundary conditions

    a = complex(zeros(Nx));
    a(1) = 2*(1 + ((m*dx^2)/h^2)*V(1) - (1i*2*m*dx^2)/(h*dt));
    for k = 2:Nx-1
        a(k) = 2*(1 + ((m*dx^2)/h^2)*V(k) - (1i*2*m*dx^2)/(h*dt)) - 1/a(k-1);
    end
    
    % Time loop
    
    Omega = zeros(Nx, Nt);
    b = complex(zeros(Nx, Nt));
    for n = 2:Nt
        for k = 2:Nx-1
            Omega(k, n) = -psi(k-1, n-1)+2*( (1i*2*m*(dx^2))/(h*dt) + 1 + (m*(dx^2)/(h^2))*V(k) ) * psi(k, n-1) - psi(k+1, n-1);
        end
    
        b(1, n) = Omega(1, n);
        for k = 2:Nx-1
            b(k, n) = b(k-1, n)/a(k-1) + Omega(k, n);
        end
    
        for k = Nx-1:-1:2
            psi(k, n) = 1/a(k) * (psi(k+1, n) - b(k, n));
        end
        
        %Normalization
        psi(:, n) = sqrt(1./trapz(psi(:, n).*conj(psi(:, n)))) * psi(:, n);
    end

    %integral_value = trapz(psi.*conj(psi));
end

function psi0 = initial_psi0(x, sigma, q, x_0, h)
    prefactor = 1/(sqrt(sigma*sqrt(pi)));
    exp1 = exp(-(x-x_0).^2/(2*sigma^2));
    exp2 = exp((1i/h) * q * x);

    psi0 = prefactor * (exp1 .* exp2);

    Norm = sqrt(1/trapz(psi0.*conj(psi0)));
    psi0 = Norm*psi0;

    % integral_value = trapz(psi0.*conj(psi0));
end

function [V1, V2] = potential(V0, x, a, b, d)
    V1 = V0 * (heaviside(x - a) - heaviside(x - (a + d)));
    V2 = V1 + V0 * (heaviside(x - b) - heaviside(x - (b + d)));
end

%function v = calculate_velocity(x, psi_square, dt)
%    max_index = find(psi_square(:) == max(psi_square(:)));
%    v = (x(max_index+1) - x(max_index-1)) / (dt);
%end

% function d = calculate_width(x, psi_square)
%     half_max = max(psi_square) / 2;
%     indices = find(psi_square >= half_max); % Find the indices where the probability density is greater or equal to half-max
%     d = x(indices(end)) - x(indices(1));
% end