%% Demo: Comprehensive Analysis of Curved Shells
% This demo performs multiple types of analysis on different shell geometries:
% - Modal analysis of multiple shell types
% - Comparison of natural frequencies
% - Dynamic response analysis
% - Stress analysis
% - Energy distribution

clear;
close all;
addpath('..');

%% Initialize Shell Structures
fprintf('=== Comprehensive Shell Analysis ===\n');
fprintf('Creating shell structures...\n');

% Common parameters
t = 0.002;      % Thickness (m)
E = 70e9;       % Young's modulus (Pa) - Aluminum
nu = 0.33;      % Poisson's ratio
rho = 2700;     % Density (kg/m³)
zeta = 0.02;    % Damping ratio

% Create cylindrical shell
cyl_params = struct('R', 0.3, 'L', 0.6, 'W', 0.4, 't', t, 'E', E, 'nu', nu, 'rho', rho, 'zeta', zeta);
cylinder = CurvedShellAnalysis.CylindricalSurface(cyl_params);

% Create spherical shell
sph_params = struct('R', 0.3, 't', t, 'phi', pi/3, 'E', E, 'nu', nu, 'rho', rho);
sphere = CurvedShellAnalysis.SphericalSurface(sph_params);

% Create conical shell
con_params = struct('R1', 0.3, 'R2', 0.1, 'L', 0.5, 't', t, 'E', E, 'nu', nu, 'rho', rho);
cone = CurvedShellAnalysis.ConicalSurface(con_params);

%% Modal Analysis
fprintf('Performing modal analysis...\n');

% Analyze each shell
shells = {cylinder, sphere, cone};
shell_names = {'Cylinder', 'Sphere', 'Cone'};
modes = cell(1,3);
freqs = cell(1,3);

figure('Position', [100 100 1200 800]);
for i = 1:3
    modal = CurvedShellAnalysis.ModalAnalysis(shells{i}, 10);  % Compute 10 modes
    modal.analyze();
    modes{i} = modal;
    freqs{i} = modal.frequencies;
    
    % Plot first 4 modes
    for j = 1:4
        subplot(3,4,(i-1)*4+j);
        modal.plotMode(j);
        title(sprintf('%s Mode %d (%.1f Hz)', shell_names{i}, j, freqs{i}(j)));
    end
end
saveas(gcf, '../docs/images/comparative_modes.png');

%% Frequency Comparison
fprintf('Comparing natural frequencies...\n');

figure('Position', [100 100 800 400]);
hold on;
markers = {'o-', 's-', 'd-'};
for i = 1:3
    plot(1:6, freqs{i}(1:6), markers{i}, 'LineWidth', 1.5);
end
grid on;
xlabel('Mode Number');
ylabel('Natural Frequency (Hz)');
title('Natural Frequency Comparison');
legend(shell_names);
saveas(gcf, '../docs/images/frequency_comparison.png');

%% Dynamic Response Analysis
fprintf('Calculating dynamic response...\n');

% Apply impulse load and calculate response
t = linspace(0, 0.1, 1000);  % Time vector
F = zeros(size(t));
F(1) = 1000;  % Impulse at t=0

figure('Position', [100 100 1200 400]);
for i = 1:3
    subplot(1,3,i);
    response = modes{i}.calculateTransientResponse(F, t);
    plot(t*1000, response);
    grid on;
    xlabel('Time (ms)');
    ylabel('Displacement (m)');
    title(sprintf('%s Response', shell_names{i}));
end
saveas(gcf, '../docs/images/dynamic_response.png');

%% Strain Energy Distribution
fprintf('Analyzing strain energy distribution...\n');

figure('Position', [100 100 1200 400]);
for i = 1:3
    subplot(1,3,i);
    energy = modes{i}.calculateModalEnergy(1);  % First mode
    shells{i}.plotScalarField(energy, sprintf('%s Strain Energy', shell_names{i}));
end
saveas(gcf, '../docs/images/strain_energy.png');

%% Statistical Analysis
fprintf('Performing statistical analysis...\n');

% Calculate MAC matrices
figure('Position', [100 100 1200 400]);
for i = 1:3
    subplot(1,3,i);
    mac = modes{i}.calculateMAC(6);  % First 6 modes
    imagesc(mac);
    colorbar;
    title(sprintf('%s MAC Matrix', shell_names{i}));
    xlabel('Mode Number');
    ylabel('Mode Number');
end
saveas(gcf, '../docs/images/mac_matrices.png');

% Print summary statistics
fprintf('\nSummary Statistics:\n');
fprintf('Shell Type\tMean Freq (Hz)\tStd Dev\tMax Freq (Hz)\n');
for i = 1:3
    f = freqs{i}(1:6);
    fprintf('%s\t\t%.1f\t\t%.1f\t%.1f\n', ...
        shell_names{i}, mean(f), std(f), max(f));
end

fprintf('\nAnalysis completed! Results saved in docs/images directory\n');
