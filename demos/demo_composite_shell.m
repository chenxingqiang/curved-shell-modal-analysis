%% Composite Shell Analysis Demo
% Author: Xingqiang Chen
% Date: 2024-12-23

% Clear workspace
clear; clc; close all;
addpath('..');

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
params.R = 0.5;     % Radius of curvature (m)
params.zeta = 0.02; % Damping ratio
params.nx = 10;     % Number of elements in x direction (reduced for stability)
params.ny = 10;     % Number of elements in y direction (reduced for stability)

%% Define Composite Properties
% Layer properties
n_layers = 4;  % Number of layers
layer_thickness = 0.002;  % Thickness per layer (m) (increased for stability)
layer_angles = [0, 45, -45, 90];  % Fiber orientations

% Material properties for each layer (carbon fiber/epoxy composite)
E1 = 140e9;    % Longitudinal modulus (Pa)
E2 = 10e9;     % Transverse modulus (Pa)
G12 = 5e9;     % Shear modulus (Pa)
nu12 = 0.3;    % Major Poisson's ratio
rho = 1600;    % Density (kg/m³)

% Create layer properties
params.Layers = cell(1, n_layers);
for i = 1:n_layers
    params.Layers{i} = struct('E1', E1, 'E2', E2, 'G12', G12, ...
                             'nu12', nu12, 'rho', rho);
end
params.LayerAngles = layer_angles;
params.LayerThickness = layer_thickness * ones(1, n_layers);

% Create composite shell
fprintf('Creating composite cylindrical shell...\n');
shell = CurvedShellAnalysis.CompositeShell(params);

%% Modal Analysis
fprintf('Performing modal analysis...\n');

% Create modal analysis object
n_modes_requested = 6;
modal = CurvedShellAnalysis.ModalAnalysis(shell, n_modes_requested);
modal.analyze();

% Plot mode shapes
n_modes = length(modal.frequencies);
n_rows = ceil(n_modes/3);
n_cols = min(3, n_modes);

figure('Name', 'Composite Shell Mode Shapes', 'Position', [100 100 1200 800]);
for i = 1:n_modes
    subplot(n_rows, n_cols, i);
    modal.plotMode(i);
    title(sprintf('Mode %d: %.2f Hz', i, full(modal.frequencies(i))));
    colorbar;
end
saveas(gcf, '../docs/images/composite_modes.png');

%% Frequency Response Analysis
fprintf('Calculating frequency response...\n');

% Define frequency range
f = logspace(0, 3, 1000);  % 1 Hz to 1 kHz
omega = 2*pi*f;

% Calculate FRF (Frequency Response Function)
H = zeros(size(f));
for i = 1:length(f)
    % Dynamic stiffness matrix
    D = -omega(i)^2 + 2*1i*shell.zeta*omega(i) + 1;
    % FRF magnitude
    H(i) = 1/abs(D);
end

% Plot FRF
figure('Name', 'Frequency Response', 'Position', [100 100 800 400]);
semilogx(f, 20*log10(H), 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Response Function');
saveas(gcf, '../docs/images/frequency_response.png');

%% MAC Analysis
fprintf('Performing MAC analysis...\n');

% Calculate MAC matrix
MAC = zeros(n_modes);
for i = 1:n_modes
    for j = 1:n_modes
        MAC(i,j) = abs(modal.modes(:,i)'*modal.modes(:,j))^2 / ...
                   ((modal.modes(:,i)'*modal.modes(:,i))*(modal.modes(:,j)'*modal.modes(:,j)));
    end
end

% Plot MAC matrix
figure('Name', 'MAC Matrix', 'Position', [100 100 600 500]);
imagesc(MAC);
colormap('jet');
colorbar;
axis equal tight;
xlabel('Mode Number');
ylabel('Mode Number');
title('Modal Assurance Criterion (MAC) Matrix');
saveas(gcf, '../docs/images/mac_matrices.png');

%% Statistical Analysis
fprintf('Performing statistical analysis...\n');

% Calculate statistics
mean_freq = mean(modal.frequencies);
std_freq = std(modal.frequencies);
max_freq = max(modal.frequencies);

% Print summary
fprintf('\nSummary Statistics:\n');
fprintf('Mean Frequency: %.1f Hz\n', full(mean_freq));
fprintf('Standard Deviation: %.1f Hz\n', full(std_freq));
fprintf('Maximum Frequency: %.1f Hz\n', full(max_freq));

% Save statistics to file
if ~exist('../ModalResults', 'dir')
    mkdir('../ModalResults');
end

fid = fopen('../ModalResults/NaturalFrequencies.txt', 'w');
fprintf(fid, 'Modal Analysis Results for Composite Shell\n');
fprintf(fid, '=====================================\n\n');
fprintf(fid, 'Shell Parameters:\n');
fprintf(fid, '- Length: %.3f m\n', shell.L);
fprintf(fid, '- Width: %.3f m\n', shell.W);
fprintf(fid, '- Radius: %.3f m\n', shell.R);
fprintf(fid, '- Total Thickness: %.3f m\n', shell.TotalThickness);
fprintf(fid, '- Number of Layers: %d\n', n_layers);
fprintf(fid, '- Layer Angles: [%s]\n', num2str(layer_angles));
fprintf(fid, '\nMaterial Properties:\n');
fprintf(fid, '- E1: %.2e Pa\n', E1);
fprintf(fid, '- E2: %.2e Pa\n', E2);
fprintf(fid, '- G12: %.2e Pa\n', G12);
fprintf(fid, '- nu12: %.2f\n', nu12);
fprintf(fid, '- Density: %.1f kg/m³\n', rho);
fprintf(fid, '\nNatural Frequencies:\n');
for i = 1:n_modes
    fprintf(fid, 'Mode %d: %.2f Hz\n', i, full(modal.frequencies(i)));
end
fprintf(fid, '\nStatistical Summary:\n');
fprintf(fid, '- Mean Frequency: %.1f Hz\n', full(mean_freq));
fprintf(fid, '- Standard Deviation: %.1f Hz\n', full(std_freq));
fprintf(fid, '- Maximum Frequency: %.1f Hz\n', full(max_freq));
fclose(fid);

fprintf('\nAnalysis completed! Results saved in docs/images directory\n');
