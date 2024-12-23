classdef ThermalAnalysis < handle
    % ThermalAnalysis Class for thermal analysis of curved shells
    
    properties
        Surface         % Surface object
        alpha1         % Longitudinal thermal expansion coefficient
        alpha2         % Transverse thermal expansion coefficient
        dT            % Temperature change
        deformation   % Thermal deformation field
    end
    
    methods
        function obj = ThermalAnalysis(surface)
            % Constructor for ThermalAnalysis
            obj.Surface = surface;
        end
        
        function setThermalProperties(obj, alpha1, alpha2)
            % Set thermal expansion coefficients
            obj.alpha1 = alpha1;
            obj.alpha2 = alpha2;
        end
        
        function setTemperatureChange(obj, dT)
            % Set temperature change
            obj.dT = dT;
        end
        
        function analyze(obj)
            % Perform thermal analysis
            fprintf('开始热分析...\n');
            
            % Get mesh from surface
            mesh = obj.Surface.getMesh();
            
            % Get material properties
            E = obj.Surface.E;
            nu = obj.Surface.nu;
            t = obj.Surface.t;
            
            % Get mesh dimensions
            [ny, nx] = size(mesh.X);
            ndof = 3;  % 3 DOF per node (u, v, w)
            
            % Initialize matrices
            n = nx * ny * ndof;
            K = sparse(n, n);
            F = zeros(n, 1);
            
            % Assembly
            fprintf('组装刚度矩阵和热力载荷...\n');
            for i = 1:nx-1
                for j = 1:ny-1
                    % Element nodes
                    nodes = [j+(i-1)*ny, j+1+(i-1)*ny, j+1+i*ny, j+i*ny];
                    
                    % Element coordinates
                    xe = mesh.X(nodes([1,2,3,4]));
                    ye = mesh.Y(nodes([1,2,3,4]));
                    ze = mesh.Z(nodes([1,2,3,4]));
                    
                    % Element matrices
                    [Ke, Fe] = obj.elementMatrices(xe, ye, ze, E, nu, t);
                    
                    % Global DOFs
                    dofs = [];
                    for node = nodes
                        dofs = [dofs, (node-1)*ndof+1:node*ndof];
                    end
                    
                    % Assembly
                    K(dofs,dofs) = K(dofs,dofs) + Ke;
                    F(dofs) = F(dofs) + Fe;
                end
            end
            
            % Apply boundary conditions
            fprintf('应用边界条件...\n');
            % Fix edges
            fixed_nodes = [1:ny, ny*(nx-1)+1:ny*nx];  % First and last rows
            fixed_dofs = [];
            for node = fixed_nodes
                fixed_dofs = [fixed_dofs, (node-1)*ndof+1:node*ndof];
            end
            free_dofs = setdiff(1:n, fixed_dofs);
            
            % Add small stiffness to prevent singularity
            K = K + sparse(1:n, 1:n, 1e-6*max(diag(K)), n, n);
            
            % Solve system
            fprintf('求解位移场...\n');
            u = zeros(n, 1);
            u(free_dofs) = K(free_dofs,free_dofs) \ F(free_dofs);
            
            % Store deformation field
            obj.deformation = reshape(u, ndof, []);
        end
        
        function [Ke, Fe] = elementMatrices(obj, xe, ye, ze, E, nu, t)
            % Calculate element matrices for thermal analysis
            
            % Element size
            dx = diff(xe([1,2]));
            dy = diff(ye([1,4]));
            
            % Number of DOFs per node
            ndof = 3;
            n = 4 * ndof;  % Total DOFs for element (4 nodes x 3 DOF)
            
            % Initialize matrices
            Ke = zeros(n, n);
            Fe = zeros(n, 1);
            
            % Gauss points
            gp = [-1/sqrt(3), 1/sqrt(3)];
            w = [1, 1];
            
            % Integration
            for i = 1:2
                for j = 1:2
                    % Shape functions
                    xi = gp(i);
                    eta = gp(j);
                    
                    N = 0.25 * [(1-xi)*(1-eta);
                               (1+xi)*(1-eta);
                               (1+xi)*(1+eta);
                               (1-xi)*(1+eta)];
                    
                    dNdxi = 0.25 * [-(1-eta), -(1-xi);
                                    (1-eta), -(1+xi);
                                    (1+eta),  (1+xi);
                                   -(1+eta),  (1-xi)];
                    
                    % Jacobian
                    xe_vec = xe(:);
                    ye_vec = ye(:);
                    J = [dNdxi(:,1)'*xe_vec  dNdxi(:,1)'*ye_vec;
                         dNdxi(:,2)'*xe_vec  dNdxi(:,2)'*ye_vec];
                    detJ = det(J);
                    invJ = inv(J);
                    
                    % B matrix (strain-displacement)
                    B = zeros(3, n);
                    for k = 1:4
                        dN = [dNdxi(k,1)*invJ(1,1) + dNdxi(k,2)*invJ(2,1);
                              dNdxi(k,1)*invJ(1,2) + dNdxi(k,2)*invJ(2,2)];
                        
                        idx = (k-1)*ndof + (1:ndof);
                        B(:,idx) = [dN(1)    0     0;
                                  0       dN(2)    0;
                                  dN(2)    dN(1)   0];
                    end
                    
                    % D matrix (constitutive)
                    D = zeros(3,3);
                    D(1:2,1:2) = E/(1-nu^2) * [1, nu; nu, 1];
                    D(3,3) = E/(2*(1+nu));
                    
                    % Thermal strain vector
                    eps_th = [obj.alpha1; obj.alpha2; 0] * obj.dT;
                    
                    % Element matrices
                    Ke = Ke + t * B' * D * B * detJ * w(i) * w(j);
                    Fe = Fe + t * B' * D * eps_th * detJ * w(i) * w(j);
                end
            end
            
            % Add bending contribution
            D_bend = E*t^3/(12*(1-nu^2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
            kappa_th = [obj.alpha1; obj.alpha2; 0] * obj.dT / t;
            for k = 1:4
                idx = (k-1)*ndof + 3;  % Only affect w-DOF
                Ke(idx,idx) = Ke(idx,idx) + D_bend(1,1);
                Fe(idx) = Fe(idx) + D_bend(1,1) * kappa_th(1);
            end
        end
        
        function plotDeformation(obj)
            % Plot thermal deformation
            mesh = obj.Surface.getMesh();
            [ny, nx] = size(mesh.X);
            
            % Extract displacements
            u = reshape(obj.deformation(1,:), ny, nx);
            v = reshape(obj.deformation(2,:), ny, nx);
            w = reshape(obj.deformation(3,:), ny, nx);
            
            % Plot deformed surface
            surf(mesh.X + u, mesh.Y + v, mesh.Z + w);
            shading interp;
            colormap('jet');
            colorbar;
            axis equal;
            view(45, 30);
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title(sprintf('Thermal Deformation (ΔT = %.1f°C)', obj.dT));
        end
    end
end
