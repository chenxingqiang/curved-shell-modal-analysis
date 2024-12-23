classdef advanced_tests < matlab.unittest.TestCase
    % Test suite for advanced analysis capabilities
    
    properties
        TestTolerance = 1e-6;
    end
    
    methods(TestMethodSetup)
        function setupTest(testCase)
            % Add package to path
            addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
        end
    end
    
    methods(Test)
        %% Multi-Scale Analysis Tests
        function testMultiScaleModel(testCase)
            % Test multi-scale model creation and basic properties
            
            % Create RVE parameters
            rve_params = struct();
            rve_params.size = [100e-6, 100e-6];
            rve_params.fiber_radius = 5e-6;
            rve_params.fiber_volume = 0.6;
            rve_params.mesh_size = 1e-6;
            
            % Material properties
            rve_params.fiber = struct('E', 230e9, 'nu', 0.2);
            rve_params.matrix = struct('E', 3.5e9, 'nu', 0.35);
            
            % Create model
            model = CurvedShellAnalysis.MultiScaleModel(rve_params);
            
            % Verify RVE creation
            testCase.verifyClass(model, 'CurvedShellAnalysis.MultiScaleModel');
            testCase.verifyEqual(model.FiberVolumeFraction, 0.6, 'RelTol', testCase.TestTolerance);
        end
        
        function testMicroDamage(testCase)
            % Test microscale damage evolution
            
            % Create damage model
            damage = CurvedShellAnalysis.MicroDamage();
            damage.enableFiberFailure('max_strain', 0.02);
            damage.enableMatrixCracking('max_stress', 50e6);
            
            % Test strain-based failure
            strain = 0.025;  % Above failure strain
            testCase.verifyTrue(damage.checkFiberFailure(strain));
            
            % Test stress-based failure
            stress = 60e6;  % Above failure stress
            testCase.verifyTrue(damage.checkMatrixFailure(stress));
        end
        
        function testHomogenization(testCase)
            % Test homogenization procedure
            
            % Create RVE
            rve_params = struct();
            rve_params.size = [100e-6, 100e-6];
            rve_params.fiber_volume = 0.6;
            
            % Material properties
            E_f = 230e9;  % Fiber modulus
            E_m = 3.5e9;  % Matrix modulus
            nu_f = 0.2;   % Fiber Poisson's ratio
            nu_m = 0.35;  % Matrix Poisson's ratio
            
            rve_params.fiber = struct('E', E_f, 'nu', nu_f);
            rve_params.matrix = struct('E', E_m, 'nu', nu_m);
            
            % Create homogenizer
            homogenizer = CurvedShellAnalysis.MicroScaleHomogenizer(rve_params);
            
            % Apply strain
            strain = [0.001; 0; 0; 0; 0; 0];
            [C_eff, ~] = homogenizer.homogenize(strain);
            
            % Verify effective properties using Voigt bound
            E_upper = E_f * 0.6 + E_m * 0.4;
            testCase.verifyLessThanOrEqual(C_eff(1,1), E_upper);
        end
        
        %% Fluid-Structure Interaction Tests
        function testFSIModel(testCase)
            % Test FSI model setup and basic properties
            
            % Create FSI parameters
            params = struct();
            params.L = [1.0, 0.5, 0.3];
            params.shell_t = 0.002;
            params.fluid_rho = 1.225;
            params.fluid_mu = 1.81e-5;
            
            % Create FSI model
            fsi = CurvedShellAnalysis.FSIModel(params);
            
            % Verify fluid properties
            testCase.verifyEqual(fsi.FluidDensity, 1.225, 'RelTol', testCase.TestTolerance);
            testCase.verifyEqual(fsi.FluidViscosity, 1.81e-5, 'RelTol', testCase.TestTolerance);
        end
        
        function testFluidSolver(testCase)
            % Test fluid solver
            
            % Create FSI model
            params = struct();
            params.L = [1.0, 0.5, 0.3];
            params.inlet_velocity = 10;
            
            fsi = CurvedShellAnalysis.FSIModel(params);
            
            % Set fluid properties
            fluid = struct();
            fluid.mesh_size = 0.01;
            fluid.time_step = 0.001;
            fsi.setFluidProperties(fluid);
            
            % Solve single time step
            [p, v] = fsi.solveFluid();
            
            % Verify mass conservation
            div_v = fsi.computeDivergence(v);
            testCase.verifyLessThan(max(abs(div_v)), 1e-10);
        end
        
        %% Topology Optimization Tests
        function testTopologyOptimization(testCase)
            % Test topology optimization setup and constraints
            
            % Create optimization parameters
            params = struct();
            params.R = 0.4;
            params.t = 0.003;
            params.vol_frac = 0.3;
            
            % Create optimization model
            topo = CurvedShellAnalysis.TopologyOptimization(params);
            
            % Set parameters
            opt_params = struct();
            opt_params.filter_radius = 0.02;
            opt_params.penalty = 3;
            topo.setParameters(opt_params);
            
            % Initialize design
            x = topo.initializeDesign('uniform');
            
            % Verify volume constraint
            vol = topo.computeVolume(x);
            testCase.verifyEqual(vol, params.vol_frac, 'RelTol', 0.01);
        end
        
        function testSensitivityAnalysis(testCase)
            % Test sensitivity analysis for topology optimization
            
            % Create simple test problem
            params = struct();
            params.L = [0.1, 0.1];
            params.vol_frac = 0.4;
            
            topo = CurvedShellAnalysis.TopologyOptimization(params);
            
            % Add load and boundary conditions
            topo.addSupport('left', 'fixed');
            topo.addPointLoad([0.1, 0.05], [0, -1000]);
            
            % Initialize design
            x = topo.initializeDesign('uniform');
            
            % Compute sensitivities
            [u, f] = topo.analyze(x);
            dc = topo.computeSensitivities(u);
            
            % Verify sensitivity properties
            testCase.verifySize(dc, size(x));
            testCase.verifyTrue(all(isfinite(dc(:))));
        end
        
        %% Dynamic Contact Tests
        function testContactDetection(testCase)
            % Test contact detection algorithm
            
            % Create contact surfaces
            params = struct();
            params.R1 = 0.3;
            params.R2 = 0.25;
            params.gap = 0.01;
            
            assembly = CurvedShellAnalysis.ShellAssembly(params);
            
            % Create contact pair
            contact = CurvedShellAnalysis.ContactPair();
            contact.setMasterSurface(assembly.Shell1);
            contact.setSlaveSurface(assembly.Shell2);
            
            % Test point projection
            point = [0.3, 0, 0];
            [dist, proj] = contact.projectPoint(point);
            
            % Verify projection properties
            testCase.verifyGreaterThanOrEqual(dist, 0);
            testCase.verifySize(proj, [1, 3]);
        end
        
        function testFrictionForces(testCase)
            % Test friction force computation
            
            % Create contact detector
            master = CurvedShellAnalysis.Surface();
            slave = CurvedShellAnalysis.Surface();
            detector = CurvedShellAnalysis.ContactDetector(master, slave);
            
            % Set friction coefficient
            detector.friction = 0.2;
            
            % Test friction force
            v_t = [0.1, 0, 0];
            f_n = [0, 0, -1000];
            f_t = detector.computeFrictionForce(v_t, f_n);
            
            % Verify Coulomb friction law
            testCase.verifyEqual(norm(f_t), 0.2*norm(f_n), 'RelTol', testCase.TestTolerance);
        end
        
        %% Time Integration Tests
        function testAdaptiveTimeIntegration(testCase)
            % Test adaptive time integration
            
            % Create test system
            M = eye(2);
            C = zeros(2);
            K = [2 -1; -1 2];
            F = [0; 1];
            
            integrator = CurvedShellAnalysis.AdaptiveTimeIntegrator(M, C, K, F);
            
            % Set parameters
            integrator.tolerance = 1e-6;
            integrator.dt_min = 1e-4;
            integrator.dt_max = 1e-2;
            
            % Solve
            tspan = [0, 0.1];
            [t, u, v, a] = integrator.solve(tspan);
            
            % Verify solution properties
            testCase.verifyGreaterThan(length(t), 1);
            testCase.verifyEqual(size(u,1), 2);
            testCase.verifyTrue(all(isfinite(u(:))));
        end
        
        function testTimeStepControl(testCase)
            % Test time step control mechanism
            
            % Create oscillator system
            M = eye(1);
            C = zeros(1);
            K = [100];  % Stiff system
            F = [1];
            
            integrator = CurvedShellAnalysis.AdaptiveTimeIntegrator(M, C, K, F);
            
            % Set parameters
            integrator.tolerance = 1e-6;
            integrator.dt_min = 1e-6;
            integrator.dt_max = 1e-3;
            
            % Take single step
            [u1, v1, a1] = integrator.step(0, 0, 0, 1e-3);
            [u2, v2, a2] = integrator.step(0, 0, 0, 5e-4);
            
            % Verify error estimation
            error = norm(u2 - u1)/norm(u2);
            testCase.verifyLessThan(error, integrator.tolerance);
        end
    end
end
