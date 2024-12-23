%% Composite Shell Analysis Demo
% Author: Xingqiang Chen
% Date: 2024-12-23

% Clear workspace
clear; clc; close all;

% Create results directories
if ~exist('results', 'dir')
    mkdir('results');
end
if ~exist('../docs/images', 'dir')
    mkdir('../docs/images');
end

%% Define Shell Parameters
fprintf('\n=== Composite Shell Analysis ===\n');

% Basic shell parameters
params = struct();
params.L = 0.4;     % Length (m)
params.W = 0.3;     % Width (m)
params.t = 0.002;   % Total thickness (m)
params.R = 0.5;     % Radius of curvature (m)
params.nu = 0.3;    % Poisson's ratio
params.zeta = 0.02; % Damping ratio

% Create cylindrical shell
fprintf('Creating composite cylindrical shell...\n');
shell = CurvedShellAnalysis.CylindricalSurface(params);

%% Define Composite Properties
% Layer properties
n_layers = 4;  % Number of layers
layer_thickness = params.t / n_layers;
layer_angles = [0, 45, -45, 90];  % Fiber orientations

% Material properties for each layer
E1 = 138e9;  % Longitudinal modulus (Pa)
E2 = 8.96e9; % Transverse modulus (Pa)
G12 = 7.1e9; % Shear modulus (Pa)
nu12 = 0.3;  % Major Poisson's ratio
rho = 1600;  % Density (kg/m³)

% Set material properties
material = struct();
material.E1 = E1;
material.E2 = E2;
material.G12 = G12;
material.nu12 = nu12;
material.rho = rho;
material.layup = layer_angles;
material.thickness = layer_thickness * ones(1, n_layers);

shell.setMaterial(material);

%% Modal Analysis
fprintf('Performing modal analysis...\n');

% Create modal analysis object
modal = CurvedShellAnalysis.ModalAnalysis(shell, 6);
modal.analyze();

% Plot mode shapes
figure('Name', 'Composite Shell Mode Shapes', 'Position', [100 100 1200 800]);
for i = 1:6
    subplot(2,3,i);
    modal.plotMode(i);
    title(sprintf('Mode %d: %.2f Hz', i, modal.frequencies(i)));
    colorbar;
end
saveas(gcf, '../docs/images/composite_modes.png');

%% Thermal Analysis
fprintf('Performing thermal analysis...\n');

% Define thermal load
dT = 50;  % Temperature change (°C)
alpha1 = -0.1e-6;  % Longitudinal thermal expansion coefficient
alpha2 = 25.6e-6;  % Transverse thermal expansion coefficient

% Calculate thermal deformation
thermal = CurvedShellAnalysis.ThermalAnalysis(shell);
thermal.setThermalProperties(alpha1, alpha2);
thermal.setTemperatureChange(dT);
thermal.analyze();

% Plot thermal deformation
figure('Name', 'Thermal Deformation', 'Position', [100 100 800 600]);
thermal.plotDeformation();
colorbar;
title(sprintf('Thermal Deformation (ΔT = %.1f°C)', dT));
saveas(gcf, '../docs/images/thermal_deformation.png');

%% Save Results
fprintf('\nAnalysis completed! Results saved in docs/images directory\n');
