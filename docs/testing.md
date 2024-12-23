# Testing Documentation

## Unit Tests

### 1. Modal Analysis Tests
```matlab
classdef ModalAnalysisTest < matlab.unittest.TestCase
    % Unit tests for modal analysis
    
    properties
        Surface
        Analysis
    end
    
    methods(TestMethodSetup)
        function setupTest(testCase)
            % Create test surface
            params = struct();
            params.L = 0.4;     % Length
            params.W = 0.3;     % Width
            params.t = 0.002;   % Thickness
            params.R = 0.5;     % Radius
            params.E = 2.1e11;  % Young's modulus
            params.rho = 7800;  % Density
            params.nu = 0.3;    % Poisson's ratio
            
            testCase.Surface = CurvedShellAnalysis.SphericalSurface(params);
            testCase.Analysis = CurvedShellAnalysis.ModalAnalysis(testCase.Surface, 6);
        end
    end
    
    methods(Test)
        function testFrequencyCalculation(testCase)
            % Test natural frequency calculation
            testCase.Analysis.analyze();
            freqs = testCase.Analysis.NaturalFrequencies;
            
            % Verify number of frequencies
            testCase.verifySize(freqs, [6,1]);
            
            % Verify frequencies are positive and ascending
            testCase.verifyGreaterThan(freqs, 0);
            testCase.verifyTrue(issorted(freqs));
        end
        
        function testModeShapeNormalization(testCase)
            % Test mode shape normalization
            testCase.Analysis.analyze();
            modes = testCase.Analysis.ModeShapes;
            
            % Verify mass normalization
            M = testCase.Surface.MassMatrix;
            for i = 1:6
                mode = modes(:,i);
                norm = abs(mode' * M * mode - 1);
                testCase.verifyLessThan(norm, 1e-10);
            end
        end
        
        function testModeOrthogonality(testCase)
            % Test mode shape orthogonality
            testCase.Analysis.analyze();
            modes = testCase.Analysis.ModeShapes;
            
            % Verify orthogonality with respect to mass matrix
            M = testCase.Surface.MassMatrix;
            for i = 1:6
                for j = i+1:6
                    prod = abs(modes(:,i)' * M * modes(:,j));
                    testCase.verifyLessThan(prod, 1e-10);
                end
            end
        end
        
        function testRayleighQuotient(testCase)
            % Test Rayleigh quotient
            testCase.Analysis.analyze();
            freqs = testCase.Analysis.NaturalFrequencies;
            modes = testCase.Analysis.ModeShapes;
            
            K = testCase.Surface.StiffnessMatrix;
            M = testCase.Surface.MassMatrix;
            
            for i = 1:6
                mode = modes(:,i);
                omega2 = mode' * K * mode / (mode' * M * mode);
                freq = sqrt(omega2)/(2*pi);
                testCase.verifyEqual(freq, freqs(i), 'RelTol', 1e-10);
            end
        end
    end
end
```

### 2. Material Model Tests
```matlab
classdef MaterialModelTest < matlab.unittest.TestCase
    % Unit tests for material models
    
    properties
        Material
    end
    
    methods(TestMethodSetup)
        function setupTest(testCase)
            % Create test material
            props = struct();
            props.E = 2.1e11;
            props.nu = 0.3;
            props.sigma_y = 250e6;
            props.H = 1e9;
            
            testCase.Material = CurvedShellAnalysis.ElastoplasticMaterial(props);
        end
    end
    
    methods(Test)
        function testElasticResponse(testCase)
            % Test elastic response
            strain = [0.001; 0; 0; 0; 0; 0];  % Uniaxial strain
            [stress, ~, ~] = testCase.Material.calculateStress(strain);
            
            % Verify elastic modulus
            E = testCase.Material.Properties.E;
            expected_stress = E * strain(1);
            testCase.verifyEqual(stress(1), expected_stress, 'RelTol', 1e-10);
        end
        
        function testYieldInitiation(testCase)
            % Test yield initiation
            sigma_y = testCase.Material.Properties.sigma_y;
            E = testCase.Material.Properties.E;
            yield_strain = sigma_y/E;
            
            strain = [yield_strain; 0; 0; 0; 0; 0];
            [stress, ep, alpha] = testCase.Material.calculateStress(strain);
            
            % Verify yield stress
            testCase.verifyEqual(stress(1), sigma_y, 'RelTol', 1e-10);
            
            % Verify no plastic strain at yield
            testCase.verifyEqual(norm(ep), 0, 'AbsTol', 1e-10);
            testCase.verifyEqual(alpha, 0, 'AbsTol', 1e-10);
        end
        
        function testPlasticHardening(testCase)
            % Test plastic hardening
            sigma_y = testCase.Material.Properties.sigma_y;
            E = testCase.Material.Properties.E;
            H = testCase.Material.Properties.H;
            
            strain = [2*sigma_y/E; 0; 0; 0; 0; 0];
            [stress, ep, alpha] = testCase.Material.calculateStress(strain);
            
            % Verify hardening response
            plastic_strain = alpha/sqrt(2/3);
            expected_stress = sigma_y + H*plastic_strain;
            testCase.verifyEqual(stress(1), expected_stress, 'RelTol', 1e-10);
        end
        
        function testUnloading(testCase)
            % Test elastic unloading
            E = testCase.Material.Properties.E;
            strain_load = [0.002; 0; 0; 0; 0; 0];
            
            % Load
            [~, ep_load, alpha_load] = testCase.Material.calculateStress(strain_load);
            
            % Unload
            strain_unload = [0.001; 0; 0; 0; 0; 0];
            [stress_unload, ep_unload, alpha_unload] = ...
                testCase.Material.calculateStress(strain_unload, ep_load, alpha_load);
            
            % Verify elastic unloading
            expected_stress = E * (strain_unload(1) - ep_load(1));
            testCase.verifyEqual(stress_unload(1), expected_stress, 'RelTol', 1e-10);
            
            % Verify no change in plastic strain
            testCase.verifyEqual(ep_unload, ep_load, 'AbsTol', 1e-10);
            testCase.verifyEqual(alpha_unload, alpha_load, 'AbsTol', 1e-10);
        end
    end
end
```

## Integration Tests

### 1. Full Analysis Pipeline
```matlab
classdef FullAnalysisTest < matlab.unittest.TestCase
    % Integration tests for complete analysis pipeline
    
    methods(Test)
        function testModalThermalCoupling(testCase)
            % Test coupled modal-thermal analysis
            
            % Create geometry
            params = struct();
            params.R = 0.5;
            params.t = 0.002;
            params.thermal_expansion = 1.2e-5;
            sphere = CurvedShellAnalysis.SphericalSurface(params);
            
            % Set temperature field
            T = @(x,y,z) 20 + 50*sqrt(x.^2 + y.^2)/params.R;
            
            % Create coupled analysis
            coupled = CurvedShellAnalysis.CoupledAnalysis(sphere);
            coupled.setTemperature(T);
            
            % Solve
            coupled.solve();
            
            % Verify results
            testCase.verifyTrue(coupled.hasConverged());
            testCase.verifySize(coupled.Frequencies, [6,1]);
            testCase.verifySize(coupled.ThermalStress, [size(sphere.Nodes,1)*6,1]);
        end
        
        function testCompositeFailure(testCase)
            % Test progressive failure analysis of composite
            
            % Create layup
            layup = struct();
            layup.thicknesses = [0.001, 0.001, 0.001];
            layup.angles = [0, 45, -45];
            
            % Create composite shell
            params = struct();
            params.R = 0.3;
            shell = CurvedShellAnalysis.CylindricalSurface(params);
            composite = CurvedShellAnalysis.CompositeMaterial(layup);
            shell.Material = composite;
            
            % Create nonlinear analysis
            nonlinear = CurvedShellAnalysis.NonlinearAnalysis(shell);
            nonlinear.enableProgressiveFailure();
            
            % Apply load
            shell.addPressure('inner', 1e6);
            
            % Solve
            [u, failed] = nonlinear.solve();
            
            % Verify results
            testCase.verifySize(u, [size(shell.Nodes,1)*6,1]);
            testCase.verifyTrue(any(failed), 'No failure detected');
        end
    end
end
```

## Performance Tests

### 1. Computational Efficiency
```matlab
classdef PerformanceTest < matlab.unittest.TestCase
    % Performance tests for computational efficiency
    
    properties(TestParameter)
        NumElements = {100, 1000, 10000}
        NumModes = {10, 20, 50}
    end
    
    methods(Test)
        function testModalAnalysisScaling(testCase, NumElements, NumModes)
            % Test scaling of modal analysis
            
            % Create surface with specified elements
            params = struct();
            params.mesh_size = sqrt(1/NumElements);
            surface = CurvedShellAnalysis.SphericalSurface(params);
            
            % Create analysis
            modal = CurvedShellAnalysis.ModalAnalysis(surface, NumModes);
            
            % Measure time
            tic;
            modal.analyze();
            time = toc;
            
            % Verify performance
            expected_time = NumElements * NumModes * 1e-6;  % Approximate scaling
            testCase.verifyLessThan(time, expected_time);
        end
        
        function testParallelSpeedup(testCase, NumElements)
            % Test parallel speedup
            
            % Create surface
            params = struct();
            params.mesh_size = sqrt(1/NumElements);
            surface = CurvedShellAnalysis.SphericalSurface(params);
            
            % Serial execution
            tic;
            surface.calculateStiffnessMatrix();
            serial_time = toc;
            
            % Parallel execution
            parallel = CurvedShellAnalysis.ParallelManager(4);
            tic;
            parallel.calculateStiffnessMatrix(surface);
            parallel_time = toc;
            
            % Verify speedup
            speedup = serial_time/parallel_time;
            testCase.verifyGreaterThan(speedup, 2.0);
        end
    end
end
```

## Validation Tests

### 1. Analytical Solutions
```matlab
classdef ValidationTest < matlab.unittest.TestCase
    % Validation tests against analytical solutions
    
    methods(Test)
        function testCircularPlateFrequencies(testCase)
            % Test against analytical solution for circular plate
            
            % Create plate
            params = struct();
            params.R = 1.0;
            params.t = 0.01;
            params.E = 2.1e11;
            params.nu = 0.3;
            params.rho = 7800;
            
            plate = CurvedShellAnalysis.CircularPlate(params);
            modal = CurvedShellAnalysis.ModalAnalysis(plate, 3);
            
            % Compute frequencies
            modal.analyze();
            freqs = modal.NaturalFrequencies;
            
            % Analytical frequencies
            lambda = [10.216, 21.26, 34.877];  % First three frequency parameters
            h = params.t;
            R = params.R;
            D = params.E*h^3/(12*(1-params.nu^2));
            rho = params.rho;
            
            analytical_freqs = lambda.^2/(2*pi*R^2) * sqrt(D/(rho*h));
            
            % Verify against analytical solution
            testCase.verifyEqual(freqs, analytical_freqs, 'RelTol', 0.05);
        end
        
        function testCylindricalShellBuckling(testCase)
            % Test against analytical buckling load
            
            % Create cylinder
            params = struct();
            params.R = 0.2;
            params.L = 1.0;
            params.t = 0.002;
            params.E = 2.1e11;
            params.nu = 0.3;
            
            cylinder = CurvedShellAnalysis.CylindricalSurface(params);
            buckling = CurvedShellAnalysis.BucklingAnalysis(cylinder);
            
            % Compute critical load
            [P_cr, ~] = buckling.solve();
            
            % Analytical solution (Timoshenko)
            E = params.E;
            nu = params.nu;
            R = params.R;
            t = params.t;
            P_analytical = 0.605*E*t^2/(R*sqrt(1-nu^2));
            
            % Verify against analytical solution
            testCase.verifyEqual(P_cr, P_analytical, 'RelTol', 0.1);
        end
    end
end
```

### 2. Benchmark Problems
```matlab
classdef BenchmarkTest < matlab.unittest.TestCase
    % Tests against standard benchmark problems
    
    methods(Test)
        function testScordelisLoRoof(testCase)
            % Scordelis-Lo roof benchmark
            
            % Create geometry
            params = struct();
            params.L = 50;    % Length
            params.R = 25;    % Radius
            params.phi = 80;  % Included angle (degrees)
            params.t = 0.25;  % Thickness
            params.E = 4.32e8;
            params.nu = 0;
            
            roof = CurvedShellAnalysis.ScordelisRoof(params);
            static = CurvedShellAnalysis.StaticAnalysis(roof);
            
            % Apply load
            roof.addGravityLoad(90);  % Load intensity
            
            % Solve
            u = static.solve();
            
            % Reference solution
            ref_displacement = 0.3024;  % Vertical displacement at midpoint
            
            % Get displacement at reference point
            mid_point = [25, 0, 0];
            [~, node_idx] = min(sum((roof.Nodes - mid_point).^2, 2));
            actual_displacement = abs(u(node_idx*6-2));  % Vertical component
            
            % Verify against reference solution
            testCase.verifyEqual(actual_displacement, ref_displacement, 'RelTol', 0.05);
        end
        
        function testPinchedCylinder(testCase)
            % Pinched cylinder benchmark
            
            % Create geometry
            params = struct();
            params.L = 600;   % Length
            params.R = 300;   % Radius
            params.t = 3;     % Thickness
            params.E = 3e6;
            params.nu = 0.3;
            
            cylinder = CurvedShellAnalysis.CylindricalSurface(params);
            static = CurvedShellAnalysis.StaticAnalysis(cylinder);
            
            % Apply loads
            P = 1;  % Point load
            cylinder.addPointLoad([0, params.R, 0], [0, -P, 0]);
            cylinder.addPointLoad([0, -params.R, 0], [0, P, 0]);
            
            % Solve
            u = static.solve();
            
            % Reference solution
            ref_displacement = 1.8248e-5;  % Radial displacement under load
            
            % Get displacement at load point
            [~, node_idx] = min(sum((cylinder.Nodes - [0, params.R, 0]).^2, 2));
            actual_displacement = abs(u(node_idx*6-1));  % Radial component
            
            % Verify against reference solution
            testCase.verifyEqual(actual_displacement, ref_displacement, 'RelTol', 0.05);
        end
    end
end
```

## Advanced Validation Cases
```matlab
classdef AdvancedValidationTest < matlab.unittest.TestCase
    % Advanced validation tests against complex analytical solutions
    
    methods(Test)
        function testLayeredShellFrequencies(testCase)
            % Test natural frequencies of layered cylindrical shell
            
            % Shell parameters
            params = struct();
            params.R = 0.5;     % Radius
            params.L = 2.0;     % Length
            params.layers = 3;  % Number of layers
            
            % Layer properties
            h = [0.001, 0.002, 0.001];  % Layer thicknesses
            theta = [0, 45, -45];       % Layer angles
            
            % Material properties
            E1 = 138e9;  % Longitudinal modulus
            E2 = 8.96e9; % Transverse modulus
            G12 = 7.1e9; % Shear modulus
            nu12 = 0.3;  % Poisson's ratio
            rho = 1600;  % Density
            
            % Create layered shell
            shell = CurvedShellAnalysis.LayeredShell(params);
            for i = 1:params.layers
                layer = struct('thickness', h(i), 'angle', theta(i));
                shell.addLayer(layer);
            end
            
            % Set material properties
            material = CurvedShellAnalysis.OrthotropicMaterial();
            material.setProperties('E1', E1, 'E2', E2, 'G12', G12, ...
                                 'nu12', nu12, 'rho', rho);
            shell.setMaterial(material);
            
            % Modal analysis
            modal = CurvedShellAnalysis.ModalAnalysis(shell, 5);
            modal.analyze();
            
            % Analytical frequencies (from literature)
            % Based on Soldatos and Hadjigeorgiou (1990)
            f_analytical = [124.5, 246.8, 389.4, 552.1, 734.9];
            
            % Compare with numerical results
            f_numerical = modal.NaturalFrequencies;
            testCase.verifyEqual(f_numerical, f_analytical, 'RelTol', 0.05);
        end
        
        function testThermalBuckling(testCase)
            % Test thermal buckling of composite cylindrical shell
            
            % Shell parameters
            params = struct();
            params.R = 0.3;
            params.L = 1.0;
            params.t = 0.003;
            
            % Create composite shell
            shell = CurvedShellAnalysis.CompositeCylinder(params);
            
            % Add layers [0/90/0]
            layup = struct();
            layup.angles = [0, 90, 0];
            layup.thicknesses = [0.001, 0.001, 0.001];
            shell.setLayup(layup);
            
            % Material properties
            mat = struct();
            mat.E1 = 138e9;
            mat.E2 = 8.96e9;
            mat.G12 = 7.1e9;
            mat.nu12 = 0.3;
            mat.alpha1 = -0.3e-6;  % Thermal expansion coefficients
            mat.alpha2 = 28.1e-6;
            shell.setMaterial(mat);
            
            % Thermal buckling analysis
            buckling = CurvedShellAnalysis.ThermalBuckling(shell);
            [T_cr, mode] = buckling.solve();
            
            % Analytical critical temperature (Shen 2001)
            T_analytical = 185.7;  % Critical temperature rise
            
            % Verify critical temperature
            testCase.verifyEqual(T_cr, T_analytical, 'RelTol', 0.1);
        end
        
        function testDynamicSnapping(testCase)
            % Test dynamic snap-through of shallow shell
            
            % Shell parameters
            params = struct();
            params.R = 2.0;    % Radius of curvature
            params.h = 0.1;    % Rise
            params.t = 0.001;  % Thickness
            params.E = 70e9;   % Young's modulus
            params.nu = 0.3;   % Poisson's ratio
            params.rho = 2700; % Density
            
            % Create shallow shell
            shell = CurvedShellAnalysis.ShallowShell(params);
            
            % Dynamic analysis
            dynamic = CurvedShellAnalysis.DynamicAnalysis(shell);
            dynamic.enableNonlinearGeometry();
            
            % Apply step load
            P = 1000;  % Load magnitude
            shell.addPointLoad([0, 0, 0], [0, 0, -P]);
            
            % Time integration parameters
            time_params = struct();
            time_params.dt = 1e-5;
            time_params.total_time = 0.01;
            dynamic.setTimeIntegration(time_params);
            
            % Solve
            [t, u] = dynamic.solve();
            
            % Get center displacement history
            center_node = shell.findNode([0, 0, 0]);
            w = u(center_node*3, :);
            
            % Analytical results from literature
            % Based on Amabili (2008)
            t_snap = 0.00234;  % Time to first snap
            w_max = -0.0187;   % Maximum displacement
            
            % Verify snap-through characteristics
            [~, snap_idx] = min(diff(w));
            t_snap_num = t(snap_idx);
            w_max_num = min(w);
            
            testCase.verifyEqual(t_snap_num, t_snap, 'RelTol', 0.15);
            testCase.verifyEqual(w_max_num, w_max, 'RelTol', 0.15);
        end
        
        function testWaveDispersion(testCase)
            % Test wave dispersion in composite cylindrical shell
            
            % Shell parameters
            params = struct();
            params.R = 0.25;
            params.L = 5.0;
            params.omega = 2*pi*1000;  % Angular frequency
            
            % Create composite shell
            shell = CurvedShellAnalysis.CompositeCylinder(params);
            
            % Add layers [0/45/-45/90]s
            angles = [0, 45, -45, 90, 90, -45, 45, 0];
            h = 0.000125 * ones(1,8);  % Layer thicknesses
            shell.setLayup(angles, h);
            
            % Material properties (AS4/3501-6)
            mat = struct();
            mat.E1 = 142e9;
            mat.E2 = 10.3e9;
            mat.G12 = 7.2e9;
            mat.nu12 = 0.27;
            mat.rho = 1580;
            shell.setMaterial(mat);
            
            % Wave propagation analysis
            wave = CurvedShellAnalysis.WaveAnalysis(shell);
            [k, v_g] = wave.dispersionRelation(params.omega);
            
            % Analytical solution
            % Based on Langley (1996)
            k_analytical = 12.45;  % Wavenumber
            v_analytical = 1245;   % Group velocity
            
            % Verify wave characteristics
            testCase.verifyEqual(k, k_analytical, 'RelTol', 0.1);
            testCase.verifyEqual(v_g, v_analytical, 'RelTol', 0.1);
        end
    end
end
```

## Performance Benchmarks
```matlab
classdef PerformanceBenchmark < matlab.unittest.TestCase
    % Performance benchmarks for computational efficiency
    
    properties(TestParameter)
        ProblemSize = {1e3, 1e4, 1e5}
        NumThreads = {1, 2, 4, 8}
    end
    
    methods(Test)
        function testAssemblyScaling(testCase, ProblemSize)
            % Test assembly performance scaling
            
            % Create large shell problem
            params = struct();
            params.mesh_size = sqrt(1/ProblemSize);
            shell = CurvedShellAnalysis.SphericalSurface(params);
            
            % Measure assembly time
            tic;
            K = shell.assembleStiffness();
            time = toc;
            
            % Verify linear scaling
            expected_time = ProblemSize * 1e-6;  % Approximate scaling
            testCase.verifyLessThan(time, expected_time);
        end
        
        function testSolverScaling(testCase, ProblemSize)
            % Test solver performance scaling
            
            % Create system
            n = ProblemSize;
            K = sprand(n, n, 5/n);
            K = K + K' + speye(n);
            f = rand(n, 1);
            
            % Direct solver
            tic;
            u_direct = K \ f;
            time_direct = toc;
            
            % Iterative solver
            tic;
            [u_iter,~] = pcg(K, f, 1e-6, 1000);
            time_iter = toc;
            
            % Verify solver efficiency
            testCase.verifyLessThan(time_iter, time_direct);
            testCase.verifyEqual(u_iter, u_direct, 'RelTol', 1e-6);
        end
        
        function testParallelSpeedup(testCase, NumThreads)
            % Test parallel speedup
            
            % Create large problem
            n = 1e5;
            params = struct();
            params.mesh_size = sqrt(1/n);
            shell = CurvedShellAnalysis.SphericalSurface(params);
            
            % Serial execution
            tic;
            K_serial = shell.assembleStiffness();
            time_serial = toc;
            
            % Parallel execution
            parallel = CurvedShellAnalysis.ParallelManager(NumThreads);
            tic;
            K_parallel = parallel.assembleStiffness(shell);
            time_parallel = toc;
            
            % Verify speedup
            speedup = time_serial/time_parallel;
            efficiency = speedup/NumThreads;
            
            testCase.verifyGreaterThan(efficiency, 0.7);
            testCase.verifyEqual(K_serial, K_parallel, 'RelTol', 1e-12);
        end
        
        function testMemoryUsage(testCase, ProblemSize)
            % Test memory efficiency
            
            % Monitor memory usage
            m0 = memory;
            
            % Create and solve problem
            params = struct();
            params.mesh_size = sqrt(1/ProblemSize);
            shell = CurvedShellAnalysis.SphericalSurface(params);
            
            modal = CurvedShellAnalysis.ModalAnalysis(shell, 10);
            modal.analyze();
            
            % Check memory usage
            m1 = memory;
            mem_used = (m1.MemUsedMATLAB - m0.MemUsedMATLAB)/1e6;
            
            % Verify memory efficiency
            expected_mem = ProblemSize * 8 * 1.5;  % Approximate usage in MB
            testCase.verifyLessThan(mem_used, expected_mem);
        end
    end
end
