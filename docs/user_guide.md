# User Guide

## Table of Contents
1. [Getting Started](#getting-started)
2. [Basic Analysis](#basic-analysis)
3. [Advanced Features](#advanced-features)
4. [Tips and Best Practices](#tips-and-best-practices)
5. [Troubleshooting](#troubleshooting)

## Getting Started

### Installation
1. System Requirements:
   - MATLAB R2020b or later
   - Parallel Computing Toolbox (optional)
   - Optimization Toolbox (optional)

2. Installation Steps:
   ```matlab
   % Add package to MATLAB path
   addpath(genpath('path/to/CurvedShellAnalysis'));
   
   % Verify installation
   help CurvedShellAnalysis
   ```

### Basic Configuration
```matlab
% Set default parameters
params = struct();
params.mesh_size = 0.01;     % Element size
params.solver_type = 'direct';  % 'direct' or 'iterative'
params.plot_results = true;   % Enable visualization

% Configure parallel computing
if license('test', 'Distrib_Computing_Toolbox')
    parpool('local');  % Start parallel pool
end
```

## Basic Analysis

### Modal Analysis Tutorial
```matlab
% 1. Define geometry parameters
params = struct();
params.L = 0.4;     % Length
params.W = 0.3;     % Width
params.t = 0.002;   % Thickness
params.R = 0.5;     % Radius

% 2. Define material properties
params.E = 2.1e11;  % Young's modulus
params.rho = 7800;  % Density
params.nu = 0.3;    % Poisson's ratio

% 3. Create surface
sphere = CurvedShellAnalysis.SphericalSurface(params);

% 4. Create modal analysis
modal = CurvedShellAnalysis.ModalAnalysis(sphere, 6);

% 5. Run analysis
modal.analyze();

% 6. Plot results
modal.plotMode(1);  % Plot first mode
```

### Static Analysis Tutorial
```matlab
% 1. Create surface (as above)
sphere = CurvedShellAnalysis.SphericalSurface(params);

% 2. Apply boundary conditions
sphere.addConstraint('bottom', 'fixed');
sphere.addLoad('top', [0, 0, -1000]);  % 1000N vertical load

% 3. Create static analysis
static = CurvedShellAnalysis.StaticAnalysis(sphere);

% 4. Solve
displacement = static.solve();

% 5. Post-process
static.plotDeformation();
static.plotStress('vonMises');
```

### Thermal Analysis Tutorial
```matlab
% 1. Define thermal properties
params.alpha = 1.2e-5;  % Thermal expansion
params.k = 45;         % Thermal conductivity
params.c = 460;        % Specific heat

% 2. Create thermal-structural analysis
thermal = CurvedShellAnalysis.ThermalStructuralAnalysis(sphere);

% 3. Set temperature field
T = @(x,y,z) 20 + 50*sqrt(x.^2 + y.^2)/params.R;
thermal.setTemperature(T);

% 4. Solve coupled problem
thermal.solve();

% 5. Visualize results
thermal.plotTemperature();
thermal.plotThermalStress();
```

## Advanced Features

### Composite Material Analysis
```matlab
% 1. Define composite layup
layup = struct();
layup.thicknesses = [0.001, 0.001, 0.001];  % Layer thicknesses
layup.angles = [0, 45, -45];                 % Fiber angles
layup.materials = {'carbon', 'carbon', 'carbon'};

% 2. Create composite material
composite = CurvedShellAnalysis.CompositeMaterial(layup);

% 3. Assign to surface
sphere.Material = composite;

% 4. Add failure criteria
composite.addFailureCriterion('tsai_wu');
composite.addFailureCriterion('puck');

% 5. Analyze
analysis = CurvedShellAnalysis.StaticAnalysis(sphere);
analysis.solve();

% 6. Check failure
composite.checkFailure();
composite.plotFailureIndices();
```

### Shape Optimization
```matlab
% 1. Define optimization problem
problem = struct();
problem.objective = 'minimizeWeight';
problem.constraints = {'maxStress', 'minFrequency'};
problem.variables = {'thickness', 'radius'};
problem.bounds = [0.001, 0.005; 0.3, 0.7];

% 2. Create optimizer
opt = CurvedShellAnalysis.AdvancedOptimization(problem, 'genetic');

% 3. Run optimization
opt.optimize();

% 4. Visualize results
opt.plotConvergence();
opt.plotPareto();  % For multi-objective
```

### Parallel Computing
```matlab
% 1. Initialize parallel manager
pm = CurvedShellAnalysis.ParallelManager(4, 100, '4GB');

% 2. Add analysis tasks
for i = 1:100
    params = generateParameters(i);
    surface = CurvedShellAnalysis.SphericalSurface(params);
    analysis = CurvedShellAnalysis.ModalAnalysis(surface, 6);
    pm.addTask('modal_analysis', analysis);
end

% 3. Execute in parallel
results = pm.executeTasks();

% 4. Process results
processResults(results);
```

## Tips and Best Practices

### Mesh Refinement
1. Start with a coarse mesh for initial analysis
2. Perform mesh convergence study
3. Refine mesh in areas of interest
4. Use adaptive meshing for efficiency

### Solver Selection
- Direct solvers: Best for small to medium problems
- Iterative solvers: Better for large problems
- GPU acceleration: Beneficial for dense matrices

### Memory Management
1. Clear large arrays when not needed
2. Use sparse matrices when possible
3. Implement out-of-core solutions for very large problems

### Result Validation
1. Compare with analytical solutions when available
2. Check energy balance
3. Verify boundary condition effects
4. Monitor convergence metrics

## Troubleshooting

### Common Issues and Solutions

#### Convergence Problems
```matlab
% 1. Adjust convergence parameters
analysis.Options.tolerance = 1e-6;
analysis.Options.max_iterations = 100;

% 2. Use line search
analysis.Options.line_search = true;

% 3. Implement load stepping
analysis.Options.num_steps = 10;
```

#### Memory Issues
```matlab
% 1. Use sparse matrices
analysis.Options.use_sparse = true;

% 2. Implement domain decomposition
analysis.Options.num_domains = 4;

% 3. Enable out-of-core solution
analysis.Options.out_of_core = true;
```

#### Numerical Instabilities
```matlab
% 1. Add numerical damping
analysis.Options.alpha = 0.1;

% 2. Use stabilization
analysis.Options.stabilization = true;

% 3. Switch to mixed formulation
analysis.Options.formulation = 'mixed';
```

### Error Messages

1. "Singular stiffness matrix":
   - Check boundary conditions
   - Verify material properties
   - Look for disconnected elements

2. "Convergence failed":
   - Increase maximum iterations
   - Reduce load step size
   - Check for material instabilities

3. "Out of memory":
   - Use sparse matrices
   - Enable parallel computing
   - Implement domain decomposition
