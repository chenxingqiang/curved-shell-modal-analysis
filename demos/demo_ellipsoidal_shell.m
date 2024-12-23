%% Demo: Modal Analysis of an Ellipsoidal Shell
% This demo shows modal analysis of an ellipsoidal shell including:
% - Natural frequencies and mode shapes
% - Random vibration analysis
% - Stress analysis

clear;
close all;
addpath('..');

%% Create an ellipsoidal shell
fprintf('=== Ellipsoidal Shell Analysis ===\n');
fprintf('Creating ellipsoidal shell...\n');

% Shell parameters
a = 0.4;        % Semi-major axis (m)
b = 0.3;        % Semi-minor axis (m)
t = 0.002;      % Thickness (m)
phi = pi/2;     % Opening angle (radians)

% Create ellipsoidal surface
params = struct('a', a, 'b', b, 't', t, 'phi', phi);
shell = CurvedShellAnalysis.EllipsoidalSurface(params);

% Material properties (Steel)
E = 210e9;      % Young's modulus (Pa)
nu = 0.3;       % Poisson's ratio
rho = 7800;     % Density (kg/m³)
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
saveas(gcf, '../docs/images/ellipsoidal_modes.png');

%% Random Vibration Analysis
fprintf('Performing random vibration analysis...\n');

% Define PSD of base excitation (white noise)
f = linspace(0, 200, 1000);
Sxx = ones(size(f)) * 1e-6;  % PSD magnitude (g²/Hz)

% Calculate response PSD
node = shell.getNumberOfNodes();
[Syy, f] = modal.calculateResponsePSD(node, f, Sxx);

% Plot PSD
figure('Position', [100 100 800 400]);
loglog(f, Sxx, 'b-', 'LineWidth', 1.5);
hold on;
loglog(f, Syy, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('PSD (g²/Hz)');
legend('Input PSD', 'Response PSD');
title('Random Vibration Analysis');
saveas(gcf, '../docs/images/ellipsoidal_psd.png');

%% Stress Analysis
fprintf('Calculating modal stresses...\n');

% Calculate modal stresses for first mode
mode_number = 1;
[sigma_xx, sigma_yy, tau_xy] = modal.calculateModalStress(mode_number);

% Plot von Mises stress
figure('Position', [100 100 800 600]);
subplot(2,2,1);
shell.plotScalarField(sigma_xx);
title('\sigma_{xx}');
colorbar;

subplot(2,2,2);
shell.plotScalarField(sigma_yy);
title('\sigma_{yy}');
colorbar;

subplot(2,2,3);
shell.plotScalarField(tau_xy);
title('\tau_{xy}');
colorbar;

% Von Mises stress
sigma_vm = sqrt(sigma_xx.^2 + sigma_yy.^2 - sigma_xx.*sigma_yy + 3*tau_xy.^2);
subplot(2,2,4);
shell.plotScalarField(sigma_vm);
title('\sigma_{VM}');
colorbar;

saveas(gcf, '../docs/images/ellipsoidal_stress.png');

fprintf('\nAnalysis completed! Results saved in docs/images directory\n');
