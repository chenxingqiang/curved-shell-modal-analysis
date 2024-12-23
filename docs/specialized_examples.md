# Specialized Examples

## Advanced Shell Configurations

### 1. Stiffened Shell Analysis

```matlab
% Create stiffened cylindrical shell
params = struct();
params.R = 0.5;      % Radius
params.L = 2.0;      % Length
params.t = 0.002;    % Shell thickness
params.E = 2.1e11;   % Young's modulus
params.nu = 0.3;     % Poisson's ratio
params.rho = 7800;   % Density

% Create base cylinder
cylinder = CurvedShellAnalysis.CylindricalSurface(params);

% Add stiffeners
stiffener_props = struct();
stiffener_props.height = 0.05;    % Stiffener height
stiffener_props.width = 0.005;    % Stiffener width
stiffener_props.spacing = 0.2;    % Longitudinal spacing

% Add longitudinal stiffeners
num_long = 8;  % Number of longitudinal stiffeners
for i = 1:num_long
    angle = 2*pi*(i-1)/num_long;
    stiffener = CurvedShellAnalysis.Stiffener('longitudinal', stiffener_props);
    stiffener.setPosition('angle', angle);
    cylinder.addStiffener(stiffener);
end

% Add ring stiffeners
ring_props = stiffener_props;
ring_props.spacing = 0.4;  % Ring spacing
cylinder.addRingStiffeners(ring_props);

% Modal analysis of stiffened shell
modal = CurvedShellAnalysis.ModalAnalysis(cylinder, 10);
modal.analyze();

% Plot mode shapes with stiffeners
modal.plotMode(1, 'ShowStiffeners', true);
```

### 2. Multi-Layer Composite Dome

```matlab
% Create composite spherical dome
params = struct();
params.R = 1.0;      % Radius
params.angle = 60;   % Dome angle
params.t = 0.01;     % Total thickness

% Define composite layup
layup = struct();
layup.materials = {'carbon', 'glass', 'carbon'};
layup.thicknesses = [0.003, 0.004, 0.003];
layup.angles = [0, 45, -45];

% Material properties
carbon = struct('E1', 230e9, 'E2', 15e9, 'G12', 27e9, 'nu12', 0.2, ...
                'Xt', 2100e6, 'Xc', 1650e6, 'Yt', 61e6, 'Yc', 130e6, 'S', 48e6);
glass = struct('E1', 80e9, 'E2', 15e9, 'G12', 17e9, 'nu12', 0.25, ...
               'Xt', 1500e6, 'Xc', 1200e6, 'Yt', 50e6, 'Yc', 110e6, 'S', 42e6);

% Create composite dome
dome = CurvedShellAnalysis.SphericalDome(params);
composite = CurvedShellAnalysis.CompositeMaterial(layup);
composite.setMaterialProperties('carbon', carbon);
composite.setMaterialProperties('glass', glass);
dome.Material = composite;

% Progressive failure analysis
nonlinear = CurvedShellAnalysis.NonlinearAnalysis(dome);
nonlinear.enableProgressiveFailure();
nonlinear.addFailureCriteria({'hashin', 'puck'});

% Apply pressure load
dome.addPressure('inner', 1e6);  % 1 MPa internal pressure

% Solve with load stepping
nonlinear.Options.num_steps = 50;
nonlinear.solve();

% Visualize results
nonlinear.plotFailureProgression();
nonlinear.plotInterlaminarStress();
nonlinear.plotDelamination();
```

### 3. Thermal-Acoustic Analysis

```matlab
% Create coupled thermal-acoustic analysis
params = struct();
params.L = [1.0, 0.8, 0.6];  % Cavity dimensions
params.t = 0.002;            % Shell thickness
params.fluid_density = 1.225;  % Air density
params.sound_speed = 343;     % Speed of sound
params.thermal_conductivity = 45;
params.specific_heat = 460;
params.thermal_expansion = 1.2e-5;

% Create coupled system
system = CurvedShellAnalysis.CoupledSystem(params);
system.addAcousticCavity(params.L);
system.addThermalBoundary('walls', 20);  % 20°C wall temperature
system.addAcousticSource([0.5, 0.4, 0.3], 1000);  % Source location and frequency

% Set coupling parameters
coupling = struct();
coupling.thermal_acoustic = 0.1;  % Coupling coefficient
coupling.acoustic_thermal = 0.2;
system.setCoupling(coupling);

% Frequency response analysis
freq_range = 0:1:1000;
response = system.frequencyResponse(freq_range);

% Plot results
system.plotAcousticPressure();
system.plotTemperatureDistribution();
system.plotCoupledResponse();
```

### 4. Shape Memory Alloy Shell

```matlab
% Create shell with SMA material
params = struct();
params.R = 0.3;      % Radius
params.L = 0.6;      % Length
params.t = 0.001;    % Thickness

% SMA material properties
sma = struct();
sma.E_austenite = 70e9;     % Austenite Young's modulus
sma.E_martensite = 30e9;    % Martensite Young's modulus
sma.nu = 0.33;              % Poisson's ratio
sma.transformation_strain = 0.05;
sma.As = 273 + 50;          % Austenite start temperature
sma.Af = 273 + 60;          % Austenite finish temperature
sma.Ms = 273 + 40;          % Martensite start temperature
sma.Mf = 273 + 30;          % Martensite finish temperature

% Create SMA material model
material = CurvedShellAnalysis.SMAMaterial(sma);

% Create shell with SMA
cylinder = CurvedShellAnalysis.CylindricalSurface(params);
cylinder.Material = material;

% Thermomechanical analysis
analysis = CurvedShellAnalysis.SMAAnalysis(cylinder);

% Define temperature cycle
time = 0:0.1:10;
temp = @(t) 273 + 30 + 40*sin(2*pi*0.1*t);
analysis.setTemperature(temp);

% Apply mechanical load
cylinder.addPressure('inner', 1e5);

% Solve
analysis.solve(time);

% Plot results
analysis.plotPhaseTransformation();
analysis.plotStrainTemperature();
analysis.plotStressStrain();
```

### 5. Functionally Graded Shell

```matlab
% Create functionally graded shell
params = struct();
params.R = 0.4;      % Radius
params.t = 0.003;    % Thickness
params.n_layers = 20;  % Number of layers for FGM discretization

% Material properties at inner surface (ceramic)
ceramic = struct('E', 380e9, 'nu', 0.24, 'rho', 3800, ...
                'k', 10.4, 'alpha', 7.4e-6, 'c', 875);

% Material properties at outer surface (metal)
metal = struct('E', 200e9, 'nu', 0.3, 'rho', 7800, ...
              'k', 15, 'alpha', 12e-6, 'c', 460);

% Create FGM material
fgm = CurvedShellAnalysis.FGMaterial(ceramic, metal);
fgm.setGradingParameter(2.0);  % Power law parameter

% Create shell
sphere = CurvedShellAnalysis.SphericalSurface(params);
sphere.Material = fgm;

% Thermal stress analysis
thermal = CurvedShellAnalysis.ThermalAnalysis(sphere);
thermal.setTemperatureGradient(100);  % 100°C temperature difference

% Solve
thermal.solve();

% Plot results
thermal.plotTemperature();
thermal.plotThermalStress();
thermal.plotMaterialDistribution();
```

### 6. Piezoelectric Shell Control

```matlab
% Create smart shell with piezoelectric patches
params = struct();
params.L = [0.5, 0.4];  % Shell dimensions
params.t = 0.002;       % Shell thickness
params.t_piezo = 0.0005;  % Piezo thickness

% Base shell material
shell_material = struct('E', 70e9, 'nu', 0.3, 'rho', 2700);

% Piezoelectric material properties
piezo = struct();
piezo.E = 63e9;          % Young's modulus
piezo.nu = 0.3;          % Poisson's ratio
piezo.rho = 7600;        % Density
piezo.d31 = -190e-12;    % Piezoelectric constant
piezo.d33 = 390e-12;
piezo.epsilon = 1.5e-8;  % Permittivity

% Create shell with piezo patches
shell = CurvedShellAnalysis.SmartShell(params);
shell.setBaseMaterial(shell_material);

% Add piezo patches
patch_locations = [0.2, 0.2; 0.3, 0.3; 0.2, 0.3; 0.3, 0.2];
for i = 1:size(patch_locations, 1)
    patch = CurvedShellAnalysis.PiezoPatch(piezo);
    patch.setLocation(patch_locations(i,:));
    shell.addPiezoPatch(patch);
end

% Create controller
controller = CurvedShellAnalysis.PiezoController(shell);
controller.setControlLaw('LQR');
controller.setTargetMode(1);  % Control first mode

% Dynamic analysis with control
dynamic = CurvedShellAnalysis.DynamicAnalysis(shell);
dynamic.enableControl(controller);
dynamic.setInitialDisturbance([0.001, 0, 0]);  % Initial displacement

% Solve
time = 0:0.001:1;
response = dynamic.solve(time);

% Plot results
dynamic.plotControlledResponse();
dynamic.plotVoltageHistory();
dynamic.plotModalAmplitudes();
```

These specialized examples demonstrate advanced applications of the shell analysis package, including complex material models, multi-physics coupling, and smart structures.
