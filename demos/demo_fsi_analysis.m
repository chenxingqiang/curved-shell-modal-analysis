%% Demo: Fluid-Structure Interaction Analysis
% This demo shows the analysis of a flexible shell in fluid flow

% Add package to path
addpath('..');

%% Create FSI Model
% Domain parameters
params = struct();
params.L = [1.0, 0.5, 0.3];  % Domain dimensions (m)
params.shell_t = 0.002;      % Shell thickness (m)
params.fluid_rho = 1.225;    % Air density (kg/m³)
params.fluid_mu = 1.81e-5;   % Air viscosity (Pa·s)
params.inlet_velocity = 10;   % Inlet velocity (m/s)

% Create FSI model
fsi = CurvedShellAnalysis.FSIModel(params);

%% Set up Fluid Domain
% Fluid mesh and solver parameters
fluid = struct();
fluid.mesh_size = 0.01;
fluid.time_step = 0.001;
fluid.turbulence_model = 'k-omega';
fsi.setFluidProperties(fluid);

%% Set up Structure
% Create flexible shell
shell = CurvedShellAnalysis.FlexibleShell(params);
shell.enableLargeDeformation();
fsi.setStructure(shell);

% Material properties (Aluminum)
material = struct();
material.E = 70e9;     % Young's modulus (Pa)
material.nu = 0.33;    % Poisson's ratio
material.rho = 2700;   % Density (kg/m³)
shell.setMaterial(material);

%% Set Coupling Parameters
coupling = struct();
coupling.relaxation = 0.7;
coupling.max_iterations = 20;
coupling.tolerance = 1e-6;
fsi.setCouplingParameters(coupling);

%% Time Integration
% Time parameters
t_start = 0;
t_end = 1.0;
dt = 0.001;
time = t_start:dt:t_end;

% Storage for results
p_hist = zeros(fsi.NumFluidNodes, length(time));
u_hist = zeros(shell.NumNodes*3, length(time));
f_hist = zeros(shell.NumNodes*3, length(time));

% Initialize figure for animation
figure('Name', 'FSI Analysis');
colormap('jet');

% Time stepping loop
for i = 1:length(time)
    % Current time
    t = time(i);
    
    % Solve coupled system
    [p, u, f] = fsi.solveTimeStep(t);
    
    % Store results
    p_hist(:,i) = p;
    u_hist(:,i) = u;
    f_hist(:,i) = f;
    
    % Plot every 10 steps
    if mod(i,10) == 0
        % Pressure field
        subplot(2,1,1);
        fsi.plotPressureField(p);
        title(sprintf('Pressure Field at t = %.3f s', t));
        colorbar;
        
        % Structure deformation
        subplot(2,1,2);
        fsi.plotDeformation(u);
        title(sprintf('Deformation at t = %.3f s', t));
        colorbar;
        
        drawnow;
    end
    
    % Check for flutter onset
    if fsi.checkFlutterOnset()
        fprintf('Flutter detected at t = %.3f s\n', t);
        break;
    end
end

%% Post-processing
% Create results directory
if ~exist('results', 'dir')
    mkdir('results');
end

% Plot pressure history at monitoring point
figure('Name', 'Pressure History');
monitor_point = fsi.getMonitorPoint();
p_monitor = p_hist(monitor_point,:);
plot(time, p_monitor);
xlabel('Time (s)');
ylabel('Pressure (Pa)');
title('Pressure History at Monitor Point');
savefig('results/pressure_history.fig');

% Plot displacement history
figure('Name', 'Displacement History');
tip_node = shell.getTipNode();
u_tip = u_hist(tip_node*3-2:tip_node*3,:);
plot(time, u_tip);
xlabel('Time (s)');
ylabel('Displacement (m)');
legend('x', 'y', 'z');
title('Tip Displacement History');
savefig('results/displacement_history.fig');

% Save results
save('results/fsi_results.mat', 'time', 'p_hist', 'u_hist', 'f_hist');

% Generate report
reportGenerator = CurvedShellAnalysis.ReportGenerator();
reportGenerator.addTimeHistory('Pressure', time, p_monitor);
reportGenerator.addTimeHistory('Displacement', time, u_tip);
reportGenerator.addContour('Final Pressure', fsi.getFluidMesh(), p);
reportGenerator.addContour('Final Deformation', shell.getMesh(), u);
reportGenerator.generate('results/fsi_report.html');
