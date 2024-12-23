classdef NonlinearAnalysis < handle
    % NONLINEARANALYSIS Class for geometric and material nonlinear analysis
    
    properties
        Surface  % Surface object
        ModalAnalysis  % Reference to modal analysis object
        LoadSteps = 10  % Number of load steps
        MaxIterations = 50  % Maximum Newton-Raphson iterations
        Tolerance = 1e-6  % Convergence tolerance
        MaterialModel  % Nonlinear material model
        OutputFolder = 'NonlinearResults'  % Output folder
    end
    
    methods
        function obj = NonlinearAnalysis(modalAnalysis)
            % Constructor
            obj.ModalAnalysis = modalAnalysis;
            obj.Surface = modalAnalysis.Surface;
            
            % Create output folder
            if ~exist(obj.OutputFolder, 'dir')
                mkdir(obj.OutputFolder);
            end
        end
        
        function setMaterialModel(obj, model_type, params)
            % Set nonlinear material model
            % Supported models: 'elastoplastic', 'hyperelastic'
            obj.MaterialModel.type = model_type;
            obj.MaterialModel.params = params;
        end
        
        function [U, stress, strain] = solveNonlinear(obj, load_vector)
            % Solve nonlinear problem using Newton-Raphson method
            fprintf('求解非线性问题...\n');
            
            % Initialize
            n = length(load_vector);
            U = zeros(n, 1);
            dU = zeros(n, 1);
            lambda = linspace(0, 1, obj.LoadSteps);
            
            % Storage for load-displacement curve
            curve.lambda = zeros(1, obj.LoadSteps);
            curve.disp = zeros(1, obj.LoadSteps);
            
            % Newton-Raphson iteration
            for step = 1:obj.LoadSteps
                fprintf('载荷步 %d/%d\n', step, obj.LoadSteps);
                F_ext = lambda(step) * load_vector;
                
                for iter = 1:obj.MaxIterations
                    % Calculate internal force and tangent stiffness
                    [F_int, K_t] = obj.calculateInternalForce(U);
                    
                    % Residual
                    R = F_ext - F_int;
                    
                    % Check convergence
                    if norm(R)/norm(F_ext) < obj.Tolerance
                        break;
                    end
                    
                    % Solve for increment
                    dU = K_t\R;
                    
                    % Update displacement
                    U = U + dU;
                    
                    fprintf('  迭代 %d: 残差 = %e\n', iter, norm(R)/norm(F_ext));
                end
                
                % Store results
                curve.lambda(step) = lambda(step);
                curve.disp(step) = max(abs(U));
            end
            
            % Calculate final stress and strain
            [stress, strain] = obj.calculateStressStrain(U);
            
            % Plot and save results
            obj.plotNonlinearResults(curve, U, stress);
        end
        
        function [F_int, K_t] = calculateInternalForce(obj, U)
            % Calculate internal force vector and tangent stiffness matrix
            [X, Y, Z] = obj.Surface.generateMesh();
            nx = obj.Surface.nx;
            ny = obj.Surface.ny;
            ndof = obj.ModalAnalysis.DOFsPerNode;
            N = nx * ny * ndof;
            
            % Initialize
            F_int = zeros(N, 1);
            K_t = sparse(N, N);
            
            % Assembly loop
            for i = 1:nx-1
                for j = 1:ny-1
                    % Get element nodes
                    nodes = [j+(i-1)*ny, j+i*ny, j+1+i*ny, j+1+(i-1)*ny];
                    
                    % Get element coordinates and displacements
                    xe = X(nodes); ye = Y(nodes); ze = Z(nodes);
                    dofs = [];
                    for n = nodes
                        dofs = [dofs, (n-1)*ndof+1:n*ndof];
                    end
                    Ue = U(dofs);
                    
                    % Calculate element contributions
                    [fe, ke] = obj.calculateElementForces(xe, ye, ze, Ue);
                    
                    % Assemble
                    F_int(dofs) = F_int(dofs) + fe;
                    K_t(dofs,dofs) = K_t(dofs,dofs) + ke;
                end
            end
        end
        
        function [fe, ke] = calculateElementForces(obj, xe, ye, ze, Ue)
            % Calculate element internal force and tangent stiffness
            % This includes both geometric and material nonlinearity
            
            % Get material properties
            E = obj.Surface.E;
            nu = obj.Surface.nu;
            
            % Calculate strain and strain increment
            [eps, deps] = obj.calculateElementStrain(xe, ye, ze, Ue);
            
            % Calculate stress based on material model
            switch obj.MaterialModel.type
                case 'elastoplastic'
                    [sigma, D] = obj.calculateElastoplasticStress(eps, deps);
                case 'hyperelastic'
                    [sigma, D] = obj.calculateHyperelasticStress(eps);
                otherwise
                    error('Unknown material model');
            end
            
            % Calculate B matrix (strain-displacement)
            B = obj.calculateBMatrix(xe, ye, ze, Ue);
            
            % Calculate element matrices
            fe = B' * sigma;
            ke = B' * D * B + obj.calculateGeometricStiffness(sigma);
        end
        
        function [sigma, D] = calculateElastoplasticStress(obj, eps, deps)
            % Calculate elastoplastic stress and material tangent
            % Using von Mises yield criterion and isotropic hardening
            
            % Material parameters
            sigma_y = obj.MaterialModel.params.yield_stress;
            H = obj.MaterialModel.params.hardening_modulus;
            E = obj.Surface.E;
            nu = obj.Surface.nu;
            
            % Elastic stiffness
            De = E/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
            
            % Trial elastic stress
            sigma_trial = De * eps;
            
            % Check yield condition
            [sigma_vm, n] = obj.calculateVonMises(sigma_trial);
            
            if sigma_vm <= sigma_y
                % Elastic response
                sigma = sigma_trial;
                D = De;
            else
                % Plastic response
                dgamma = (sigma_vm - sigma_y)/(2*G + H);
                sigma = sigma_trial - 2*G*dgamma*n;
                
                % Consistent tangent
                D = De - (4*G^2/(2*G + H)) * (n * n');
            end
        end
        
        function [sigma, D] = calculateHyperelasticStress(obj, eps)
            % Calculate hyperelastic stress and material tangent
            % Using Neo-Hookean model
            
            % Material parameters
            mu = obj.MaterialModel.params.shear_modulus;
            lambda = obj.MaterialModel.params.bulk_modulus;
            
            % Deformation gradient
            F = reshape(eps, 3, 3) + eye(3);
            C = F' * F;
            I1 = trace(C);
            J = det(F);
            
            % Second Piola-Kirchhoff stress
            S = mu*(eye(3) - C^(-1)) + lambda*log(J)*C^(-1);
            
            % Cauchy stress
            sigma = (1/J) * F * S * F';
            
            % Material tangent
            D = obj.calculateHyperelasticTangent(F, mu, lambda);
        end
        
        function [stress, strain] = calculateStressStrain(obj, U)
            % Calculate stress and strain fields from displacement
            [X, Y, Z] = obj.Surface.generateMesh();
            nx = obj.Surface.nx;
            ny = obj.Surface.ny;
            
            % Initialize fields
            stress = struct('xx', zeros(ny,nx), 'yy', zeros(ny,nx), ...
                          'xy', zeros(ny,nx), 'vm', zeros(ny,nx));
            strain = struct('xx', zeros(ny,nx), 'yy', zeros(ny,nx), ...
                          'xy', zeros(ny,nx));
            
            % Calculate element-wise quantities
            for i = 1:nx-1
                for j = 1:ny-1
                    % Element nodes and displacements
                    nodes = [j+(i-1)*ny, j+i*ny, j+1+i*ny, j+1+(i-1)*ny];
                    dofs = [];
                    for n = nodes
                        dofs = [dofs, (n-1)*5+1:n*5];
                    end
                    Ue = U(dofs);
                    
                    % Calculate strains and stresses
                    eps = obj.calculateElementStrain(X(nodes), Y(nodes), Z(nodes), Ue);
                    
                    switch obj.MaterialModel.type
                        case 'elastoplastic'
                            [sig, ~] = obj.calculateElastoplasticStress(eps, eps);
                        case 'hyperelastic'
                            [sig, ~] = obj.calculateHyperelasticStress(eps);
                    end
                    
                    % Store results
                    for k = 1:4
                        stress.xx(nodes(k)) = sig(1);
                        stress.yy(nodes(k)) = sig(2);
                        stress.xy(nodes(k)) = sig(3);
                        stress.vm(nodes(k)) = sqrt(sig(1)^2 + sig(2)^2 - sig(1)*sig(2) + 3*sig(3)^2);
                        
                        strain.xx(nodes(k)) = eps(1);
                        strain.yy(nodes(k)) = eps(2);
                        strain.xy(nodes(k)) = eps(3);
                    end
                end
            end
        end
        
        function plotNonlinearResults(obj, curve, U, stress)
            % Plot nonlinear analysis results
            
            % Load-displacement curve
            figure('Position', [100 100 800 600]);
            plot(curve.disp, curve.lambda, 'b-o', 'LineWidth', 2);
            grid on;
            title('载荷-位移曲线');
            xlabel('最大位移 (m)');
            ylabel('载荷因子');
            saveas(gcf, fullfile(obj.OutputFolder, 'LoadDisplacement.png'));
            
            % Deformed shape
            [X, Y, Z] = obj.Surface.generateMesh();
            figure('Position', [100 100 800 600]);
            
            % Original shape
            subplot(2,2,1);
            surf(X, Y, Z);
            title('原始形状');
            axis equal;
            
            % Deformed shape
            subplot(2,2,2);
            Zdef = Z + reshape(U(3:5:end), size(Z));
            surf(X, Y, Zdef);
            title('变形后形状');
            axis equal;
            
            % von Mises stress
            subplot(2,2,3);
            surf(X, Y, stress.vm);
            colormap('jet');
            colorbar;
            title('von Mises应力');
            axis equal;
            
            % Save figure
            saveas(gcf, fullfile(obj.OutputFolder, 'NonlinearResults.png'));
            
            % Export to VTK
            obj.exportToVTK(X, Y, Z, U, stress);
        end
        
        function exportToVTK(obj, X, Y, Z, U, stress)
            % Export results to VTK format
            filename = fullfile(obj.OutputFolder, 'nonlinear_results.vtk');
            
            % Open file
            fid = fopen(filename, 'w');
            
            % Write header
            fprintf(fid, '# vtk DataFile Version 3.0\n');
            fprintf(fid, 'Nonlinear Analysis Results\n');
            fprintf(fid, 'ASCII\n');
            fprintf(fid, 'DATASET STRUCTURED_GRID\n');
            fprintf(fid, 'DIMENSIONS %d %d 1\n', size(X,2), size(X,1));
            fprintf(fid, 'POINTS %d float\n', numel(X));
            
            % Write coordinates
            for i = 1:numel(X)
                fprintf(fid, '%f %f %f\n', X(i), Y(i), Z(i));
            end
            
            % Write displacements
            fprintf(fid, '\nPOINT_DATA %d\n', numel(X));
            fprintf(fid, 'VECTORS displacement float\n');
            for i = 1:numel(X)
                fprintf(fid, '%f %f %f\n', U(5*i-4), U(5*i-3), U(5*i-2));
            end
            
            % Write stresses
            fprintf(fid, '\nSCALARS vonMises float\n');
            fprintf(fid, 'LOOKUP_TABLE default\n');
            fprintf(fid, '%f\n', stress.vm(:));
            
            fclose(fid);
        end
    end
end
