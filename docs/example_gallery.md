# Example Gallery

## Modal Analysis Examples

### Spherical Shell Modes

![Spherical Shell Modes](images/sphere_modes.png)

```matlab
% Create spherical shell
params = struct('R', 0.5, 't', 0.002, 'E', 2.1e11, 'rho', 7800, 'nu', 0.3);
sphere = CurvedShellAnalysis.SphericalSurface(params);

% Modal analysis
modal = CurvedShellAnalysis.ModalAnalysis(sphere, 6);
modal.analyze();

% Plot first 6 modes
modal.plotModes(1:6, 'subplot');
```

### Cylindrical Shell Breathing Modes

![Cylinder Breathing Modes](images/cylinder_breathing.png)

```matlab
% Create cylindrical shell
params = struct('R', 0.2, 'L', 1.0, 't', 0.001);
cylinder = CurvedShellAnalysis.CylindricalSurface(params);

% Find breathing modes
modal = CurvedShellAnalysis.ModalAnalysis(cylinder, 10);
modal.findBreathingModes();
```

## Nonlinear Analysis Examples

### Large Deformation

![Large Deformation](images/large_deform.png)

```matlab
% Nonlinear analysis of cylinder under pressure
cylinder = CurvedShellAnalysis.CylindricalSurface(params);
cylinder.addPressure('inner', 1e6);  % 1 MPa internal pressure

nonlinear = CurvedShellAnalysis.NonlinearAnalysis(cylinder);
nonlinear.solve();
nonlinear.plotDeformation('scale', 2);
```

### Snap-Through Buckling

![Snap-Through](images/snap_through.png)

```matlab
% Shallow spherical cap under point load
cap = CurvedShellAnalysis.SphericalCap(params);
cap.addPointLoad('center', [0, 0, -1000]);

% Arc-length method
nonlinear = CurvedShellAnalysis.NonlinearAnalysis(cap);
nonlinear.Options.method = 'arc_length';
nonlinear.solve();
nonlinear.plotLoadDisplacement();
```

## Thermal-Structural Examples

### Thermal Stress

![Thermal Stress](images/thermal_stress.png)

```matlab
% Temperature gradient in sphere
sphere = CurvedShellAnalysis.SphericalSurface(params);
T = @(x,y,z) 20 + 50*sqrt(x.^2 + y.^2)/params.R;

thermal = CurvedShellAnalysis.ThermalStructuralAnalysis(sphere);
thermal.setTemperature(T);
thermal.solve();
thermal.plotThermalStress('vonMises');
```

### Transient Thermal

![Transient Thermal](images/transient_thermal.png)

```matlab
% Time-dependent heating
thermal.setTimeDependent(true);
thermal.setHeatSource(@(t) 1000*sin(2*pi*t));
thermal.solveTransient([0 10], 100);
thermal.animateTemperature();
```

## Composite Examples

### Layup Analysis

![Composite Layup](images/composite_layup.png)

```matlab
% Define composite layup
layup = struct();
layup.thicknesses = [0.001, 0.001, 0.001];
layup.angles = [0, 45, -45];
layup.materials = {'carbon', 'carbon', 'carbon'};

composite = CurvedShellAnalysis.CompositeMaterial(layup);
sphere.Material = composite;

% Analyze and visualize
analysis = CurvedShellAnalysis.StaticAnalysis(sphere);
analysis.solve();
composite.plotLayup();
composite.plotStressDistribution();
```

### Failure Analysis

![Failure Indices](images/failure_indices.png)

```matlab
% Add multiple failure criteria
composite.addFailureCriterion('tsai_wu');
composite.addFailureCriterion('puck');
composite.addFailureCriterion('larc');

% Check failure
composite.checkFailure();
composite.plotFailureIndices('subplot');
```

## Optimization Examples

### Shape Optimization

![Shape Optimization](images/shape_opt.png)

```matlab
% Optimize sphere radius and thickness
problem = struct();
problem.objective = 'minimizeWeight';
problem.constraints = {'maxStress', 'minFrequency'};
problem.variables = {'thickness', 'radius'};

opt = CurvedShellAnalysis.AdvancedOptimization(problem, 'genetic');
opt.optimize();
opt.plotConvergence();
```

### Multi-Objective Optimization

![Pareto Front](images/pareto_front.png)

```matlab
% Weight vs. stiffness optimization
problem.objectives = {'weight', 'compliance'};
opt = CurvedShellAnalysis.AdvancedOptimization(problem, 'nsga2');
opt.optimize();
opt.plotPareto();
```

## Advanced Visualization

### Mode Animation

![Mode Animation](images/mode_animation.gif)

```matlab
% Animate mode shape
modal = CurvedShellAnalysis.ModalAnalysis(sphere, 1);
modal.analyze();
modal.animateMode(1, 'cycles', 2, 'fps', 30);
```

### Stress Contours

![Stress Contours](images/stress_contours.png)

```matlab
% Plot stress with custom options
analysis.plotStress('vonMises', ...
    'Colormap', 'jet', ...
    'ShowMesh', true, ...
    'Deformed', true);
```

### Interactive 3D View

![Interactive View](images/interactive_view.png)

```matlab
% Create interactive visualization
viewer = CurvedShellAnalysis.ResultViewer(analysis);
viewer.addResult('displacement');
viewer.addResult('stress', 'vonMises');
viewer.show();
```
