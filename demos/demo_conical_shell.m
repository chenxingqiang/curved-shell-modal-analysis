%% Demo: Modal Analysis of a Conical Shell
% This demo shows modal analysis of a conical shell including:
% - Natural frequencies and mode shapes
% - Buckling analysis
% - Fluid-structure interaction

clear;
close all;
addpath('..');

%% Create a conical shell
fprintf('=== Conical Shell Analysis ===\n');
fprintf('Creating conical shell...\n');

% Shell parameters
R1 = 0.3;       % Base radius (m)
R2 = 0.1;       % Top radius (m)
L = 0.5;        % Length (m)
t = 0.002;      % Thickness (m)

% Create conical surface
params = struct('R1', R1, 'R2', R2, 'L', L, 't', t);
shell = CurvedShellAnalysis.ConicalSurface(params);

% Material properties (Titanium)
E = 110e9;      % Young's modulus (Pa)
nu = 0.34;      % Poisson's ratio
rho = 4500;     % Density (kg/m³)
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
saveas(gcf, '../docs/images/conical_modes.png');

%% Buckling Analysis
fprintf('Performing buckling analysis...\n');

% Create buckling analysis object
buckling = CurvedShellAnalysis.BucklingAnalysis(shell);

% Apply axial compression load
P = -1000;  % Axial force (N)
buckling.setLoad('axial', P);
buckling.analyze();

% Plot buckling modes
figure('Position', [100 100 1200 400]);
for i = 1:3
    subplot(1, 3, i);
    buckling.plotMode(i);
    title(sprintf('Buckling Mode %d (λ = %.2f)', i, buckling.eigenvalues(i)));
end
saveas(gcf, '../docs/images/conical_buckling.png');

%% Fluid-Structure Interaction
fprintf('Performing FSI analysis...\n');

% Create FSI analysis object
fsi = CurvedShellAnalysis.FSIAnalysis(shell);

% Set fluid properties (Air)
rho_f = 1.225;      % Fluid density (kg/m³)
c = 343;            % Speed of sound (m/s)
fsi.setFluidProperties(rho_f, c);

% Calculate coupled frequencies
fsi.analyze();

% Plot pressure field for first mode
figure('Position', [100 100 800 600]);
subplot(2,1,1);
fsi.plotMode(1);
title('First FSI Mode Shape');

subplot(2,1,2);
fsi.plotPressureField(1);
title('Acoustic Pressure Field');
colorbar;

saveas(gcf, '../docs/images/conical_fsi.png');

% Compare frequencies
fprintf('\nNatural Frequencies Comparison:\n');
fprintf('Mode\tIn Vacuo (Hz)\tIn Air (Hz)\tFreq. Shift\n');
for i = 1:5
    f_vac = modal.frequencies(i);
    f_fsi = fsi.frequencies(i);
    shift = (f_fsi - f_vac)/f_vac * 100;
    fprintf('%d\t%.2f\t\t%.2f\t\t%.1f%%\n', i, f_vac, f_fsi, shift);
end

fprintf('\nAnalysis completed! Results saved in docs/images directory\n');
