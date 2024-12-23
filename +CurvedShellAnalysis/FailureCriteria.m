classdef FailureCriteria < handle
    % FAILURECRITERIA Class for comprehensive failure analysis
    
    properties
        Type  % Failure criterion type
        Properties  % Material properties
        Options  % Analysis options
        Results  % Analysis results
    end
    
    methods
        function obj = FailureCriteria(type, props)
            % Constructor
            obj.Type = type;
            obj.Properties = props;
            obj.initializeOptions();
        end
        
        function initializeOptions(obj)
            % Initialize analysis options based on criterion type
            switch obj.Type
                case 'tsai_wu'
                    obj.Options.interaction = 1.0;  % Interaction term
                    
                case 'puck'
                    obj.Options.friction_angle = 30;  % Friction angle in degrees
                    obj.Options.fracture_angle = 53;  % Fracture angle in degrees
                    
                case 'larc'
                    obj.Options.shear_nonlinearity = true;
                    obj.Options.thermal_effects = true;
                    
                case 'multicontinuum'
                    obj.Options.constituent_models = {'fiber', 'matrix'};
                    obj.Options.homogenization = 'mori_tanaka';
            end
        end
        
        function [f, mode] = evaluateFailure(obj, stress, strain)
            % Evaluate failure based on selected criterion
            switch obj.Type
                case 'tsai_wu'
                    [f, mode] = obj.tsaiWuCriterion(stress);
                    
                case 'puck'
                    [f, mode] = obj.puckCriterion(stress);
                    
                case 'larc'
                    [f, mode] = obj.larcCriterion(stress, strain);
                    
                case 'multicontinuum'
                    [f, mode] = obj.multicontinuumFailure(stress, strain);
                    
                case 'cuntze'
                    [f, mode] = obj.cuntzeCriterion(stress);
                    
                case 'hashin_rotem'
                    [f, mode] = obj.hashinRotemCriterion(stress);
            end
            
            % Store results
            obj.Results.failure_index = f;
            obj.Results.failure_mode = mode;
        end
        
        function [f, mode] = tsaiWuCriterion(obj, stress)
            % Tsai-Wu failure criterion
            
            % Get material properties
            Xt = obj.Properties.Xt;  % Tensile strength in fiber direction
            Xc = obj.Properties.Xc;  % Compressive strength in fiber direction
            Yt = obj.Properties.Yt;  % Tensile strength in matrix direction
            Yc = obj.Properties.Yc;  % Compressive strength in matrix direction
            S = obj.Properties.S;    % Shear strength
            F12 = obj.Options.interaction * -1/(sqrt(Xt*Xc*Yt*Yc));  % Interaction term
            
            % Calculate strength parameters
            F1 = 1/Xt - 1/Xc;
            F2 = 1/Yt - 1/Yc;
            F11 = 1/(Xt*Xc);
            F22 = 1/(Yt*Yc);
            F66 = 1/S^2;
            
            % Calculate failure index
            f = F1*stress(1) + F2*stress(2) + ...
                F11*stress(1)^2 + F22*stress(2)^2 + F66*stress(3)^2 + ...
                2*F12*stress(1)*stress(2);
            
            % Determine failure mode
            if f >= 1
                if abs(stress(1)/Xt) > abs(stress(2)/Yt)
                    mode = 'fiber';
                else
                    mode = 'matrix';
                end
            else
                mode = 'none';
            end
        end
        
        function [f, mode] = puckCriterion(obj, stress)
            % Puck failure criterion
            
            % Get material properties
            Xt = obj.Properties.Xt;
            Xc = obj.Properties.Xc;
            Yt = obj.Properties.Yt;
            Yc = obj.Properties.Yc;
            S = obj.Properties.S;
            
            % Friction parameters
            phi = obj.Options.friction_angle * pi/180;  % Convert to radians
            theta_fp = obj.Options.fracture_angle * pi/180;
            
            % Fiber failure
            if stress(1) >= 0
                ff1 = stress(1)/Xt;
            else
                ff1 = abs(stress(1))/Xc;
            end
            
            % Matrix failure
            sigma_n = stress(2)*cos(theta_fp)^2 + ...
                     stress(3)*sin(2*theta_fp);
            tau_nt = -stress(2)*sin(2*theta_fp)/2 + ...
                     stress(3)*cos(2*theta_fp);
            tau_nl = 0;  % For 2D case
            
            if sigma_n >= 0
                % Mode A (tension)
                fm = sqrt((tau_nt/(S-sigma_n*tan(phi)))^2 + ...
                         (tau_nl/S)^2) + ...
                         sigma_n/(Yt/cos(2*theta_fp));
            else
                % Mode B/C (compression)
                fm = sqrt((tau_nt/(S-sigma_n*tan(phi)))^2 + ...
                         (tau_nl/S)^2) + ...
                         sigma_n/(Yc/cos(2*theta_fp));
            end
            
            % Get maximum failure index
            f = max(ff1, fm);
            
            % Determine failure mode
            if f >= 1
                if ff1 > fm
                    mode = 'fiber';
                else
                    if sigma_n >= 0
                        mode = 'matrix_tension';
                    else
                        mode = 'matrix_compression';
                    end
                end
            else
                mode = 'none';
            end
        end
        
        function [f, mode] = larcCriterion(obj, stress, strain)
            % LaRC failure criterion
            
            % Get material properties
            E1 = obj.Properties.E1;
            E2 = obj.Properties.E2;
            G12 = obj.Properties.G12;
            nu12 = obj.Properties.nu12;
            alpha1 = obj.Properties.alpha1;
            alpha2 = obj.Properties.alpha2;
            
            % Temperature change if thermal effects are considered
            if obj.Options.thermal_effects
                dT = obj.Properties.temperature - obj.Properties.reference_temp;
                eps_th = [alpha1*dT; alpha2*dT; 0];
            else
                eps_th = zeros(3,1);
            end
            
            % Effective strain (including thermal effects)
            eps_eff = strain - eps_th;
            
            % Calculate in-situ strengths
            Xt_is = obj.Properties.Xt;
            Xc_is = obj.Properties.Xc;
            Yt_is = obj.Properties.Yt * 1.2;  % In-situ effect
            Yc_is = obj.Properties.Yc * 1.2;
            S_is = obj.Properties.S * 1.4;
            
            % Fiber failure in tension
            if stress(1) >= 0
                FF1 = stress(1)/Xt_is;
            else
                % Fiber failure in compression (kinking)
                % Misalignment angle
                phi0 = 1.5 * pi/180;  % Initial misalignment
                
                % Rotated stresses
                sigma_m = stress(1)*cos(phi0)^2 + ...
                         stress(2)*sin(phi0)^2 + ...
                         2*stress(3)*sin(phi0)*cos(phi0);
                         
                FF1 = abs(sigma_m)/Xc_is;
            end
            
            % Matrix failure
            if stress(2) >= 0
                % Matrix tension
                if obj.Options.shear_nonlinearity
                    % Nonlinear shear behavior
                    gamma12 = strain(3);
                    tau12 = G12 * gamma12 * (1 - exp(-abs(gamma12)));
                else
                    tau12 = stress(3);
                end
                
                FM = sqrt((stress(2)/Yt_is)^2 + (tau12/S_is)^2);
            else
                % Matrix compression
                alpha0 = 53 * pi/180;  % Fracture angle
                
                % Rotated stresses
                sigma_n = stress(2)*cos(alpha0)^2 + ...
                         2*stress(3)*sin(alpha0)*cos(alpha0);
                tau_t = -stress(2)*sin(alpha0)*cos(alpha0) + ...
                        stress(3)*(cos(alpha0)^2 - sin(alpha0)^2);
                
                % Friction coefficient
                eta = 0.3;
                
                if sigma_n >= 0
                    tau_eff = sqrt(tau_t^2);
                else
                    tau_eff = sqrt(tau_t^2 + (eta*sigma_n)^2);
                end
                
                FM = tau_eff/S_is;
            end
            
            % Get maximum failure index
            f = max(FF1, FM);
            
            % Determine failure mode
            if f >= 1
                if FF1 > FM
                    if stress(1) >= 0
                        mode = 'fiber_tension';
                    else
                        mode = 'fiber_compression';
                    end
                else
                    if stress(2) >= 0
                        mode = 'matrix_tension';
                    else
                        mode = 'matrix_compression';
                    end
                end
            else
                mode = 'none';
            end
        end
        
        function [f, mode] = multicontinuumFailure(obj, stress, strain)
            % Multi-continuum failure theory
            
            % Get constituent properties
            Vf = obj.Properties.fiber_volume_fraction;
            Vm = 1 - Vf;
            
            % Perform stress decomposition
            switch obj.Options.homogenization
                case 'mori_tanaka'
                    [stress_f, stress_m] = obj.moriTanakaDecomposition(stress, Vf);
                otherwise
                    error('Unsupported homogenization scheme');
            end
            
            % Evaluate fiber failure
            ff = obj.evaluateConstituentFailure('fiber', stress_f);
            
            % Evaluate matrix failure
            fm = obj.evaluateConstituentFailure('matrix', stress_m);
            
            % Combine failure indices
            f = max(ff, fm);
            
            % Determine failure mode
            if f >= 1
                if ff > fm
                    mode = 'fiber';
                else
                    mode = 'matrix';
                end
            else
                mode = 'none';
            end
        end
        
        function [f, mode] = cuntzeCriterion(obj, stress)
            % Cuntze failure criterion (5 failure modes)
            
            % Get material properties
            Xt = obj.Properties.Xt;
            Xc = obj.Properties.Xc;
            Yt = obj.Properties.Yt;
            Yc = obj.Properties.Yc;
            S = obj.Properties.S;
            
            % Mode 1: Fiber tension
            if stress(1) >= 0
                f1 = stress(1)/Xt;
            else
                f1 = 0;
            end
            
            % Mode 2: Fiber compression
            if stress(1) < 0
                f2 = abs(stress(1))/Xc;
            else
                f2 = 0;
            end
            
            % Mode 3: Matrix tension
            if stress(2) >= 0
                f3 = sqrt((stress(2)/Yt)^2 + (stress(3)/S)^2);
            else
                f3 = 0;
            end
            
            % Mode 4: Matrix compression
            if stress(2) < 0
                f4 = sqrt((stress(2)/Yc)^2 + (stress(3)/S)^2);
            else
                f4 = 0;
            end
            
            % Mode 5: Matrix shear
            f5 = abs(stress(3))/S;
            
            % Get maximum failure index and mode
            [f, idx] = max([f1 f2 f3 f4 f5]);
            
            % Determine failure mode
            if f >= 1
                modes = {'fiber_tension', 'fiber_compression', ...
                        'matrix_tension', 'matrix_compression', 'matrix_shear'};
                mode = modes{idx};
            else
                mode = 'none';
            end
        end
        
        function [f, mode] = hashinRotemCriterion(obj, stress)
            % Hashin-Rotem failure criterion
            
            % Get material properties
            Xt = obj.Properties.Xt;
            Xc = obj.Properties.Xc;
            Yt = obj.Properties.Yt;
            Yc = obj.Properties.Yc;
            S = obj.Properties.S;
            
            % Fiber failure
            if stress(1) >= 0
                ff = stress(1)/Xt;
            else
                ff = abs(stress(1))/Xc;
            end
            
            % Matrix failure
            if stress(2) >= 0
                fm = sqrt((stress(2)/Yt)^2 + (stress(3)/S)^2);
            else
                fm = sqrt((stress(2)/Yc)^2 + (stress(3)/S)^2);
            end
            
            % Get maximum failure index
            f = max(ff, fm);
            
            % Determine failure mode
            if f >= 1
                if ff > fm
                    if stress(1) >= 0
                        mode = 'fiber_tension';
                    else
                        mode = 'fiber_compression';
                    end
                else
                    if stress(2) >= 0
                        mode = 'matrix_tension';
                    else
                        mode = 'matrix_compression';
                    end
                end
            else
                mode = 'none';
            end
        end
        
        function [stress_f, stress_m] = moriTanakaDecomposition(obj, stress, Vf)
            % Mori-Tanaka stress decomposition
            
            % Get constituent properties
            Ef = obj.Properties.E_fiber;
            Em = obj.Properties.E_matrix;
            nuf = obj.Properties.nu_fiber;
            num = obj.Properties.nu_matrix;
            
            % Calculate Eshelby tensor
            S = obj.calculateEshelbyTensor(nuf);
            
            % Calculate concentration tensors
            A = inv(eye(6) + S*((Ef/Em)-eye(6)));
            B = Vf*A + (1-Vf)*eye(6);
            
            % Calculate stresses in constituents
            stress_f = A/B * stress;
            stress_m = (eye(6) - Vf*A)/((1-Vf)*B) * stress;
        end
        
        function S = calculateEshelbyTensor(obj, nu)
            % Calculate Eshelby tensor for cylindrical fiber
            
            % Components for transversely isotropic material
            S = zeros(6,6);
            
            % Fill components (simplified for cylindrical fiber)
            S(1,1) = 5-4*nu;
            S(2,2) = S(1,1);
            S(3,3) = 0;
            S(4,4) = 1;
            S(5,5) = S(4,4);
            S(6,6) = 2*(1-nu);
            
            S = S/(8*(1-nu));
        end
        
        function f = evaluateConstituentFailure(obj, constituent, stress)
            % Evaluate failure in individual constituents
            
            switch constituent
                case 'fiber'
                    % Von Mises criterion for fiber
                    sigma_vm = sqrt(stress(1)^2 - stress(1)*stress(2) + ...
                                 stress(2)^2 + 3*stress(3)^2);
                    f = sigma_vm/obj.Properties.fiber_strength;
                    
                case 'matrix'
                    % Modified von Mises for matrix
                    sigma_vm = sqrt(stress(1)^2 - stress(1)*stress(2) + ...
                                 stress(2)^2 + 3*stress(3)^2);
                    sigma_h = (stress(1) + stress(2))/3;
                    
                    % Pressure-dependent yielding
                    f = (sigma_vm + 0.2*sigma_h)/obj.Properties.matrix_strength;
            end
        end
    end
end
