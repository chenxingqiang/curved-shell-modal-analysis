%% Demo: Topology Optimization of Curved Shell
% This demo shows topology optimization of a curved shell structure
% for minimum compliance with volume and frequency constraints

% Add package to path
addpath('..');

%% Create Optimization Model
% Shell parameters
params = struct();
params.R = 0.4;      % Radius (m)
params.t = 0.003;    % Thickness (m)
params.vol_frac = 0.3;  % Target volume fraction

% Create optimization model
topo = CurvedShellAnalysis.TopologyOptimization(params);

%% Set Optimization Parameters
opt_params = struct();
opt_params.filter_radius = 0.02;
opt_params.penalty = 3;
opt_params.max_iter = 100;
topo.setParameters(opt_params);

%% Define Problem
% Add supports
topo.addSupport('edges', 'fixed');

% Add loads
topo.addPointLoad([0.2, 0, 0], [0, -1000, 0]);
topo.addPointLoad([-0.2, 0, 0], [0, -1000, 0]);

% Set objective
topo.setObjective('compliance');

% Add constraints
topo.addConstraint('volume', params.vol_frac);
topo.addConstraint('frequency', 100, 'min');

%% Initialize Design
x = topo.initializeDesign('uniform');

% Create figure for optimization history
figure('Name', 'Optimization History');
colormap(gray);

%% Optimization Loop
history.compliance = zeros(opt_params.max_iter, 1);
history.volume = zeros(opt_params.max_iter, 1);
history.frequency = zeros(opt_params.max_iter, 1);

for iter = 1:opt_params.max_iter
    % FE analysis
    [u, f] = topo.analyze(x);
    
    % Compute objective and constraints
    compliance = topo.computeCompliance(u, f);
    volume = topo.computeVolume(x);
    freq = topo.computeFrequency(x);
    
    % Store history
    history.compliance(iter) = compliance;
    history.volume(iter) = volume;
    history.frequency(iter) = freq;
    
    % Plot current design
    subplot(2,2,1);
    topo.plotDesign(x);
    title(sprintf('Iteration %d', iter));
    
    subplot(2,2,2);
    plot(1:iter, history.compliance(1:iter));
    title('Compliance History');
    xlabel('Iteration');
    ylabel('Compliance');
    
    subplot(2,2,3);
    plot(1:iter, history.volume(1:iter));
    hold on;
    plot([1 iter], [params.vol_frac params.vol_frac], 'r--');
    hold off;
    title('Volume Fraction History');
    xlabel('Iteration');
    ylabel('Volume Fraction');
    
    subplot(2,2,4);
    plot(1:iter, history.frequency(1:iter));
    hold on;
    plot([1 iter], [100 100], 'r--');
    hold off;
    title('First Frequency History');
    xlabel('Iteration');
    ylabel('Frequency (Hz)');
    
    drawnow;
    
    % Compute sensitivities
    dc = topo.computeSensitivities(u);
    
    % Update design variables
    x_new = topo.updateDesign(x, dc);
    
    % Apply filter
    x_new = topo.filterDensities(x_new);
    
    % Check convergence
    if topo.hasConverged(x, x_new)
        break;
    end
    
    x = x_new;
end

%% Post-processing
% Create results directory
if ~exist('results', 'dir')
    mkdir('results');
end

% Final analysis
[u_final, f_final] = topo.analyze(x);

% Plot final results
figure('Name', 'Final Design');

subplot(2,2,1);
topo.plotDesign(x);
title('Optimal Design');
colorbar;

subplot(2,2,2);
topo.plotDeformation(u_final);
title('Deformation');
colorbar;

subplot(2,2,3);
topo.plotVonMises(u_final);
title('von Mises Stress');
colorbar;

subplot(2,2,4);
topo.plotModeShape(1);
title('First Mode Shape');
colorbar;

% Save results
save('results/topology_results.mat', 'x', 'history', 'u_final', 'f_final');

% Export optimal design
topo.exportSTL('results/optimal_shell.stl');

% Generate report
reportGenerator = CurvedShellAnalysis.ReportGenerator();
reportGenerator.addConvergenceHistory('Compliance', history.compliance);
reportGenerator.addConvergenceHistory('Volume', history.volume);
reportGenerator.addConvergenceHistory('Frequency', history.frequency);
reportGenerator.addContour('Final Design', topo.getMesh(), x);
reportGenerator.addContour('Deformation', topo.getMesh(), u_final);
reportGenerator.generate('results/topology_report.html');
