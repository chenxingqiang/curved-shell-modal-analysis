# Advanced Examples

## Advanced Modal Analysis

### 1. Complex Shell Intersections

```matlab
% Create intersecting cylinders
params = struct();
params.R1 = 0.3;  % First cylinder radius
params.R2 = 0.2;  % Second cylinder radius
params.L = 1.0;   % Length
params.t = 0.002; % Thickness
params.E = 2.1e11;
params.rho = 7800;
params.nu = 0.3;

% Create geometry
cyl1 = CurvedShellAnalysis.CylindricalSurface(params);
cyl2 = CurvedShellAnalysis.CylindricalSurface(params);
cyl2.rotate('y', 90);  % Rotate second cylinder

% Merge surfaces
intersection = CurvedShellAnalysis.SurfaceMerger([cyl1, cyl2]);
intersection.merge();

% Modal analysis
modal = CurvedShellAnalysis.ModalAnalysis(intersection, 10);
modal.analyze();

% Plot complex modes
modal.plotComplexModes(1:6);
```

### 2. Damped Modal Analysis

```matlab
% Define damping properties
damping = struct();
damping.alpha = 0.1;  % Mass proportional damping
damping.beta = 1e-5;  % Stiffness proportional damping
damping.modal = [0.01, 0.02];  % Modal damping ratios

% Create analysis
modal = CurvedShellAnalysis.ModalAnalysis(sphere, 10);
modal.setDamping(damping);

% Complex eigenvalue analysis
[freqs, modes, damping_ratios] = modal.analyzeComplex();

% Plot frequency response
modal.plotFRF('point', [0.1, 0, 0], ...
              'frequency_range', [0, 1000], ...
              'resolution', 1000);
```

## Advanced Nonlinear Analysis

### 1. Contact Analysis

```matlab
% Create two shells
shell1 = CurvedShellAnalysis.SphericalSurface(params1);
shell2 = CurvedShellAnalysis.SphericalSurface(params2);

% Define contact parameters
contact = struct();
contact.type = 'frictionless';
contact.normal_penalty = 1e6;
contact.detection_method = 'node_to_surface';
contact.tolerance = 1e-6;

% Create contact pair
pair = CurvedShellAnalysis.ContactPair(shell1, shell2, contact);

% Nonlinear analysis with contact
nonlinear = CurvedShellAnalysis.NonlinearAnalysis([shell1, shell2]);
nonlinear.addContactPair(pair);
nonlinear.solve();

% Visualize contact pressure
nonlinear.plotContactPressure();
nonlinear.plotContactStatus();
```

### 2. Material Nonlinearity with Damage

```matlab
% Define elastoplastic material with damage
material = struct();
material.E = 2.1e11;
material.nu = 0.3;
material.yield_stress = 250e6;
material.hardening_modulus = 1e9;
material.damage_initiation = 0.01;
material.damage_evolution = 'exponential';
material.damage_parameter = 0.1;

% Create material model
mat = CurvedShellAnalysis.AdvancedMaterial('elastoplastic_damage', material);

% Assign to surface
sphere.Material = mat;

% Nonlinear analysis
nonlinear = CurvedShellAnalysis.NonlinearAnalysis(sphere);
nonlinear.solve();

% Plot results
nonlinear.plotDamage();
nonlinear.plotPlasticStrain();
```

## Advanced Optimization

### 1. Topology Optimization

```matlab
% Define topology optimization problem
problem = struct();
problem.objective = 'minimizeCompliance';
problem.volume_fraction = 0.3;
problem.filter_radius = 0.02;
problem.penalty = 3;

% Create optimizer
topology = CurvedShellAnalysis.TopologyOptimization(sphere, problem);

% Set boundary conditions
topology.addConstraint('bottom', 'fixed');
topology.addLoad('top', [0, 0, -1000]);

% Optimize
topology.optimize('method', 'simp', ...
                 'max_iterations', 100, ...
                 'tolerance', 1e-4);

% Plot results
topology.plotDensity();
topology.plotConvergence();
```

### 2. Multi-Objective Shape Optimization

```matlab
% Define objectives
objectives = struct();
objectives.weight = @(x) calculateWeight(x);
objectives.frequency = @(x) -calculateFirstFrequency(x);
objectives.stress = @(x) calculateMaxStress(x);

% Define constraints
constraints = struct();
constraints.stress_limit = @(x) calculateMaxStress(x) - 250e6;
constraints.displacement_limit = @(x) calculateMaxDisplacement(x) - 0.001;

% Create multi-objective optimizer
opt = CurvedShellAnalysis.AdvancedOptimization('nsga2');
opt.setObjectives(objectives);
opt.setConstraints(constraints);

% Set optimization parameters
opt.Options.population_size = 100;
opt.Options.generations = 50;
opt.Options.crossover_fraction = 0.8;
opt.Options.mutation_rate = 0.1;

% Run optimization
opt.optimize();

% Analyze results
opt.plotPareto3D();
opt.plotObjectiveSpace();
opt.analyzeSensitivity();
```

## Advanced Composite Analysis

### 1. Progressive Failure Analysis

```matlab
% Define composite layup
layup = struct();
layup.thicknesses = [0.001, 0.001, 0.001, 0.001];
layup.angles = [0, 45, -45, 90];
layup.materials = {'carbon', 'carbon', 'carbon', 'carbon'};

% Create composite with multiple failure criteria
composite = CurvedShellAnalysis.CompositeMaterial(layup);
composite.addFailureCriterion('hashin');
composite.addFailureCriterion('puck');
composite.addFailureCriterion('larc');

% Progressive damage parameters
damage = struct();
damage.degradation_method = 'gradual';
damage.min_stiffness = 0.01;
damage.damage_evolution = 'exponential';

composite.setDamageModel(damage);

% Nonlinear analysis with progressive failure
nonlinear = CurvedShellAnalysis.NonlinearAnalysis(sphere);
nonlinear.Material = composite;
nonlinear.solve();

% Visualize failure progression
nonlinear.plotFailureProgression();
nonlinear.plotLoadDisplacement();
nonlinear.plotDamageContours();
```

### 2. Delamination Analysis

```matlab
% Define cohesive parameters
cohesive = struct();
cohesive.normal_strength = 30e6;
cohesive.shear_strength = 50e6;
cohesive.normal_energy = 100;
cohesive.shear_energy = 300;
cohesive.penalty_stiffness = 1e12;

% Create layup with cohesive zones
layup = CurvedShellAnalysis.CompositeMaterial(layup);
layup.addCohesiveZone(2, cohesive);  % Add between layers 2 and 3

% Nonlinear analysis
nonlinear = CurvedShellAnalysis.NonlinearAnalysis(sphere);
nonlinear.Material = layup;
nonlinear.solve();

% Plot delamination results
nonlinear.plotDelamination();
nonlinear.plotCohesiveStress();
nonlinear.plotDamageVariable();
```

## Advanced Thermal Analysis

### 1. Coupled Thermo-Mechanical Analysis

```matlab
% Define thermal properties
thermal = struct();
thermal.conductivity = 45;
thermal.specific_heat = 460;
thermal.expansion = 1.2e-5;
thermal.convection = 25;
thermal.emissivity = 0.8;

% Create coupled analysis
coupled = CurvedShellAnalysis.ThermalStructuralAnalysis(sphere);
coupled.setThermalProperties(thermal);

% Define time-dependent boundary conditions
t = 0:0.1:10;
T = @(t) 100 * sin(2*pi*0.1*t);
coupled.setTemperatureBoundary('top', T);
coupled.setConvectionBoundary('outer', 25, 20);
coupled.setRadiationBoundary('outer', 0.8, 20);

% Solve coupled problem
coupled.solveTransient(t);

% Visualize results
coupled.plotTemperatureHistory();
coupled.plotThermalStressHistory();
coupled.animateCoupledResponse();
```

### 2. Thermal Buckling Analysis

```matlab
% Setup thermal buckling analysis
buckling = CurvedShellAnalysis.ThermalBucklingAnalysis(cylinder);

% Define temperature field
T = @(x,y,z) 20 + 100*x/params.L;  % Linear temperature gradient
buckling.setTemperatureField(T);

% Solve for critical temperature
[T_cr, modes] = buckling.solve();

% Plot results
buckling.plotCriticalModes(1:3);
buckling.plotTemperatureDeformation();
```

## Advanced Specialized Examples

### 7. Multi-Scale Analysis

```matlab
% Multi-scale analysis of composite shell
params = struct();
params.R = 0.3;      % Radius
params.L = 0.6;      % Length
params.t = 0.002;    % Total thickness

% Microscale RVE properties
rve = struct();
rve.fiber_radius = 5e-6;     % Fiber radius
rve.fiber_volume = 0.6;      % Fiber volume fraction
rve.fiber_E = 230e9;         % Fiber Young's modulus
rve.fiber_nu = 0.2;          % Fiber Poisson's ratio
rve.matrix_E = 3.5e9;        % Matrix Young's modulus
rve.matrix_nu = 0.35;        % Matrix Poisson's ratio

% Create multi-scale model
model = CurvedShellAnalysis.MultiScaleModel(params);
model.setRVE(rve);

% Add damage model at microscale
damage = CurvedShellAnalysis.MicroDamage();
damage.enableFiberFailure('max_strain', 0.02);
damage.enableMatrixCracking('max_stress', 50e6);
damage.enableDelamination('critical_energy', 250);
model.addDamageModel(damage);

% Create nonlinear analysis
nonlinear = CurvedShellAnalysis.NonlinearAnalysis(model);
nonlinear.enableMultiScaleComputation();

% Apply loading
model.addPressure('inner', 2e6);

% Solve with microscale homogenization
[u, micro_damage] = nonlinear.solve();

% Visualize results
model.plotMicrostructure();
model.plotDamageEvolution();
model.plotEffectiveProperties();
```

### 8. Fluid-Structure Interaction

```matlab
% FSI analysis of shell in flow
params = struct();
params.L = [1.0, 0.5, 0.3];  % Domain dimensions
params.shell_t = 0.002;      % Shell thickness
params.fluid_rho = 1.225;    % Fluid density
params.fluid_mu = 1.81e-5;   % Fluid viscosity
params.inlet_velocity = 10;   % Inlet velocity

% Create coupled FSI model
fsi = CurvedShellAnalysis.FSIModel(params);

% Set up fluid domain
fluid = struct();
fluid.mesh_size = 0.01;
fluid.time_step = 0.001;
fluid.turbulence_model = 'k-omega';
fsi.setFluidProperties(fluid);

% Set up shell structure
shell = CurvedShellAnalysis.FlexibleShell(params);
shell.enableLargeDeformation();
fsi.setStructure(shell);

% Set coupling parameters
coupling = struct();
coupling.relaxation = 0.7;
coupling.max_iterations = 20;
coupling.tolerance = 1e-6;
fsi.setCouplingParameters(coupling);

% Time integration
time = 0:0.001:1.0;
for t = time
    % Solve coupled system
    [p, u, f] = fsi.solveTimeStep(t);
    
    % Update mesh
    fsi.updateMesh();
    
    % Check convergence
    if fsi.checkFlutterOnset()
        break;
    end
end

% Post-process results
fsi.plotPressureField();
fsi.plotDeformation();
fsi.plotVonMisesStress();
```

### 9. Topology Optimization

```matlab
% Topology optimization of curved shell
params = struct();
params.R = 0.4;      % Radius
params.t = 0.003;    % Thickness
params.vol_frac = 0.3;  % Target volume fraction

% Create optimization model
topo = CurvedShellAnalysis.TopologyOptimization(params);

% Set optimization parameters
opt_params = struct();
opt_params.filter_radius = 0.02;
opt_params.penalty = 3;
opt_params.max_iter = 100;
topo.setParameters(opt_params);

% Define loads and constraints
topo.addSupport('edges', 'fixed');
topo.addPointLoad([0.2, 0, 0], [0, -1000, 0]);

% Objective function
topo.setObjective('compliance');

% Add constraints
topo.addConstraint('volume', params.vol_frac);
topo.addConstraint('frequency', 100, 'min');

% Initialize design variables
x = topo.initializeDesign('uniform');

% Optimization loop
for iter = 1:opt_params.max_iter
    % FE analysis
    [u, f] = topo.analyze(x);
    
    % Sensitivity analysis
    dc = topo.computeSensitivities(u);
    
    % Update design variables
    x = topo.updateDesign(x, dc);
    
    % Apply filter
    x = topo.filterDensities(x);
    
    % Check convergence
    if topo.hasConverged()
        break;
    end
end

% Post-process results
topo.plotOptimalDesign();
topo.plotConvergenceHistory();
topo.exportSTL('optimal_shell.stl');
```

### 10. Dynamic Contact Analysis

```matlab
% Dynamic contact analysis of shell assembly
params = struct();
params.R1 = 0.3;     % Radius of first shell
params.R2 = 0.25;    % Radius of second shell
params.gap = 0.01;   % Initial gap
params.friction = 0.2;  % Friction coefficient

% Create shell assembly
assembly = CurvedShellAnalysis.ShellAssembly(params);

% Add contact surfaces
contact = CurvedShellAnalysis.ContactPair();
contact.setMasterSurface(assembly.Shell1);
contact.setSlaveSurface(assembly.Shell2);
contact.enableFriction(params.friction);
assembly.addContactPair(contact);

% Set up dynamic analysis
dynamic = CurvedShellAnalysis.DynamicAnalysis(assembly);
dynamic.enableNonlinearGeometry();
dynamic.enableContactDetection();

% Time integration parameters
time_params = struct();
time_params.dt = 0.0001;
time_params.total_time = 0.01;
time_params.beta = 0.25;
time_params.gamma = 0.5;
dynamic.setTimeIntegration(time_params);

% Initial conditions
v0 = zeros(size(assembly.Nodes));
v0(assembly.Shell1.NodeIds, 2) = 1.0;  % Initial velocity
dynamic.setInitialConditions([], v0);

% Solve
[u, v, a] = dynamic.solve();

% Post-process
dynamic.plotContactForces();
dynamic.plotContactStatus();
dynamic.plotDeformedShape();
```

These examples demonstrate advanced usage of the package's capabilities. Each example includes detailed code and explanations of the key features being demonstrated.
