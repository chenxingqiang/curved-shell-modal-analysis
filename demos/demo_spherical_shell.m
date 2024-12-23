%% Demo: Modal Analysis of a Spherical Shell
% This demo shows modal analysis of a spherical shell including:
% - Natural frequencies and mode shapes
% - Frequency response analysis
% - Thermal effects

clear;
close all;
addpath('..');

%% Create a spherical shell
fprintf('=== Spherical Shell Analysis ===\n');
fprintf('Creating spherical shell...\n');

% Shell parameters
R = 0.5;        % Radius (m)
t = 0.002;      % Thickness (m)
phi = pi/3;     % Opening angle (radians)

% Create spherical surface
params = struct('R', R, 't', t, 'phi', phi);
shell = CurvedShellAnalysis.SphericalSurface(params);

% Material properties (Aluminum)
E = 70e9;       % Young's modulus (Pa)
nu = 0.33;      % Poisson's ratio
rho = 2700;     % Density (kg/m³)
shell.setMaterial(E, nu, rho);

%% Modal Analysis
fprintf('Performing modal analysis...\n');
modal = CurvedShellAnalysis.ModalAnalysis(shell);
modal.analyze();

% Plot and save mode shapes
if ~exist('../docs/images', 'dir')
    mkdir('../docs/images');
end

figure('Position', [100 100 1200 800]);
for i = 1:6
    subplot(2, 3, i);
    modal.plotMode(i);
    title(sprintf('Mode %d (%.2f Hz)', i, modal.frequencies(i)));
end
saveas(gcf, '../docs/images/spherical_modes.png');

%% Frequency Response Analysis
fprintf('Calculating frequency response...\n');

% Define excitation point and response point
excitation_node = 1;
response_node = shell.getNumberOfNodes();

% Calculate FRF
f = linspace(0, 200, 1000);
[H, f] = modal.calculateFRF(excitation_node, response_node, f);

% Plot FRF
figure('Position', [100 100 800 400]);
semilogy(f, abs(H));
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (m/N)');
title('Frequency Response Function');
saveas(gcf, '../docs/images/spherical_frf.png');

%% Thermal Analysis
fprintf('Performing thermal analysis...\n');

% Create thermal analysis object
thermal = CurvedShellAnalysis.ThermalAnalysis(shell);

% Set thermal properties
alpha1 = 23e-6;  % Thermal expansion coefficient (1/K)
alpha2 = 23e-6;  % Thermal expansion coefficient (1/K)
dT = 50;        % Temperature change (K)

thermal.setThermalProperties(alpha1, alpha2);
thermal.setTemperatureChange(dT);
thermal.analyze();

% Plot thermal deformation
figure('Position', [100 100 800 600]);
thermal.plotDeformation();
saveas(gcf, '../docs/images/spherical_thermal.png');

fprintf('\nAnalysis completed! Results saved in docs/images directory\n');
