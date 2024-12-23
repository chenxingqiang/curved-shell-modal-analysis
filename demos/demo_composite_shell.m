%% Demo: Composite Shell Analysis
% This demo shows the analysis of a composite cylindrical shell under
% combined thermal and mechanical loading

% Add package to path
addpath('..');

%% Create composite shell
% Shell parameters
params = struct();
params.R = 0.3;      % Radius (m)
params.L = 0.6;      % Length (m)
params.t = 0.002;    % Total thickness (m)

% Create shell
shell = CurvedShellAnalysis.CylindricalSurface(params);

% Define composite layup [0/45/-45/90]s
angles = [0, 45, -45, 90, 90, -45, 45, 0];
thicknesses = ones(1,8) * params.t/8;

layup = struct('angles', angles, 'thicknesses', thicknesses);
shell.setLayup(layup);

%% Set material properties (Carbon/Epoxy)
material = struct();
material.E1 = 138e9;    % Longitudinal modulus (Pa)
material.E2 = 8.96e9;   % Transverse modulus (Pa)
material.G12 = 7.1e9;   % Shear modulus (Pa)
material.nu12 = 0.3;    % Major Poisson's ratio
material.rho = 1600;    % Density (kg/m³)
material.alpha1 = -0.3e-6;  % Longitudinal thermal expansion (/°C)
material.alpha2 = 28.1e-6;  % Transverse thermal expansion (/°C)

shell.setMaterial(material);

%% Modal Analysis
% Create modal analysis
modal = CurvedShellAnalysis.ModalAnalysis(shell, 10);

% Compute natural frequencies and mode shapes
modal.analyze();

% Plot first 4 mode shapes
figure('Name', 'Mode Shapes');
for i = 1:4
    subplot(2,2,i);
    modal.plotMode(i);
    title(sprintf('Mode %d: %.2f Hz', i, modal.NaturalFrequencies(i)));
end

%% Thermal-Mechanical Analysis
% Create thermal-mechanical analysis
thermal = CurvedShellAnalysis.ThermalAnalysis(shell);

% Apply temperature field
T = @(x,y,z) 20 + 50*sqrt(x.^2 + y.^2)/params.R;  % Temperature variation
thermal.setTemperature(T);

% Apply mechanical load
shell.addPressure('inner', 1e5);  % 1 bar internal pressure

% Solve
[u, stress] = thermal.solve();

% Plot results
figure('Name', 'Thermal-Mechanical Analysis');

subplot(2,1,1);
thermal.plotDeformation(u);
title('Deformed Shape');
colorbar;

subplot(2,1,2);
thermal.plotVonMises(stress);
title('von Mises Stress');
colorbar;

%% Progressive Failure Analysis
% Create nonlinear analysis
nonlinear = CurvedShellAnalysis.NonlinearAnalysis(shell);
nonlinear.enableProgressiveFailure();

% Add failure criteria
nonlinear.addFailureCriteria({'hashin', 'puck'});

% Apply load steps
load_steps = linspace(0, 2e5, 20);  % Up to 2 bar pressure
[u_hist, damage] = nonlinear.solve(load_steps);

% Plot damage evolution
figure('Name', 'Progressive Failure');
nonlinear.plotDamageEvolution(damage);
title('Damage Evolution');
xlabel('Load Step');
ylabel('Damage Variable');

% Save results
if ~exist('results', 'dir')
    mkdir('results');
end
save('results/composite_shell_results.mat', 'modal', 'thermal', 'nonlinear');

% Generate report
reportGenerator = CurvedShellAnalysis.ReportGenerator();
reportGenerator.addResults('Modal Analysis', modal);
reportGenerator.addResults('Thermal-Mechanical', thermal);
reportGenerator.addResults('Progressive Failure', nonlinear);
reportGenerator.generate('results/composite_shell_report.html');
