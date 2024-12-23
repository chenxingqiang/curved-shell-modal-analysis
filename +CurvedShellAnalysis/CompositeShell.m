classdef CompositeShell < CurvedShellAnalysis.Surface
    % CompositeShell Class for layered composite shells
    %   This class implements a layered composite shell with multiple material layers
    %   Each layer can have different material properties and orientations
    
    properties
        Layers          % Array of layer properties
        LayerAngles    % Fiber orientation angles for each layer
        LayerThickness % Thickness of each layer
        TotalThickness % Total shell thickness
        ABD            % Composite stiffness matrix (6x6)
        zeta           % Damping ratio
        R              % Radius of curvature
    end
    
    methods
        function obj = CompositeShell(params)
            % Constructor for CompositeShell
            % params: Structure with the following fields
            %   - Layers: Cell array of material properties for each layer
            %   - LayerAngles: Array of fiber orientation angles (degrees)
            %   - LayerThickness: Array of layer thicknesses
            %   - Geometry parameters (inherited from Surface)
            
            % Call superclass constructor
            obj@CurvedShellAnalysis.Surface(params);
            
            % Initialize composite properties
            if isfield(params, 'Layers')
                obj.Layers = params.Layers;
            end
            if isfield(params, 'LayerAngles')
                obj.LayerAngles = params.LayerAngles;
            end
            if isfield(params, 'LayerThickness')
                obj.LayerThickness = params.LayerThickness;
                obj.TotalThickness = sum(obj.LayerThickness);
                obj.t = obj.TotalThickness;  % Update total thickness
            end
            if isfield(params, 'zeta')
                obj.zeta = params.zeta;
            else
                obj.zeta = 0.02;  % Default damping ratio
            end
            if isfield(params, 'R')
                obj.R = params.R;
            else
                obj.R = 0.5;  % Default radius of curvature
            end
            
            % Generate mesh
            obj.generateMesh();
            
            % Calculate composite stiffness matrix
            obj.calculateABDMatrix();
        end
        
        function generateMesh(obj)
            % Generate cylindrical shell mesh
            
            % Get mesh dimensions
            nx = obj.nx + 1;  % Number of nodes in x direction
            ny = obj.ny + 1;  % Number of nodes in y direction
            
            % Create grid points
            x = linspace(0, obj.L, nx);
            y = linspace(0, obj.W, ny);
            [X, Y] = meshgrid(x, y);
            
            % Calculate theta angle for each point
            theta = (Y - obj.W/2)/obj.R;  % Angle from centerline
            
            % Calculate Z coordinates
            Z = obj.R*(1 - cos(theta));
            
            % Store mesh
            obj.mesh.X = X;
            obj.mesh.Y = obj.R*sin(theta);  % Update Y coordinates for curved surface
            obj.mesh.Z = Z;
        end
        
        function calculateABDMatrix(obj)
            % Calculate the ABD stiffness matrix for the composite layup
            n = length(obj.Layers);
            A = zeros(3,3);  % Extensional stiffness
            B = zeros(3,3);  % Coupling stiffness
            D = zeros(3,3);  % Bending stiffness
            
            % Z-coordinates of layer interfaces
            z = -obj.TotalThickness/2;
            z_coords = zeros(1,n+1);
            for i = 1:n
                z_coords(i) = z;
                z = z + obj.LayerThickness(i);
            end
            z_coords(end) = z;
            
            % Calculate stiffness contributions from each layer
            for k = 1:n
                % Get layer properties
                E1 = obj.Layers{k}.E1;
                E2 = obj.Layers{k}.E2;
                G12 = obj.Layers{k}.G12;
                nu12 = obj.Layers{k}.nu12;
                nu21 = nu12*E2/E1;  % Calculate minor Poisson's ratio
                
                % Check for physical consistency
                if nu12*nu21 >= 1
                    warning('Layer %d: Poisson''s ratios may lead to instability', k);
                    % Adjust nu12 to ensure stability
                    nu12 = 0.9/sqrt(E1/E2);
                    nu21 = nu12*E2/E1;
                end
                
                % Calculate plane stress reduced stiffness matrix
                Q = zeros(3,3);
                denom = 1 - nu12*nu21;
                Q(1,1) = E1/denom;
                Q(2,2) = E2/denom;
                Q(1,2) = nu12*E2/denom;
                Q(2,1) = Q(1,2);
                Q(3,3) = G12;
                
                % Transform to laminate coordinates
                theta = obj.LayerAngles(k)*pi/180;  % Convert to radians
                c = cos(theta);
                s = sin(theta);
                
                T = [c^2, s^2, 2*c*s;
                     s^2, c^2, -2*c*s;
                     -c*s, c*s, c^2-s^2];
                
                Qbar = T'*Q*T;
                
                % Add layer contribution to ABD matrix
                z1 = z_coords(k);
                z2 = z_coords(k+1);
                
                % Extensional stiffness (A)
                A = A + Qbar*(z2-z1);
                
                % Coupling stiffness (B)
                B = B + Qbar*(z2^2-z1^2)/2;
                
                % Bending stiffness (D)
                D = D + Qbar*(z2^3-z1^3)/3;
            end
            
            % Assemble ABD matrix
            obj.ABD = [A, B; B, D];
            
            % Add small regularization for numerical stability
            eps = 1e-6*trace(obj.ABD)/size(obj.ABD,1);
            obj.ABD = obj.ABD + eps*eye(size(obj.ABD));
            
            % Check for positive definiteness
            try
                chol(obj.ABD);
            catch
                warning('ABD matrix is not positive definite. Adding regularization...');
                % Add stronger regularization
                eps = 1e-3*trace(obj.ABD)/size(obj.ABD,1);
                obj.ABD = obj.ABD + eps*eye(size(obj.ABD));
            end
        end
        
        function Q = getReducedStiffness(obj, E1, E2, G12, nu12)
            % Calculate reduced stiffness matrix for a single layer
            nu21 = nu12*E2/E1;
            Q = zeros(3,3);
            Q(1,1) = E1/(1-nu12*nu21);
            Q(1,2) = nu12*E2/(1-nu12*nu21);
            Q(2,1) = Q(1,2);
            Q(2,2) = E2/(1-nu12*nu21);
            Q(3,3) = G12;
        end
        
        function T = getTransformationMatrix(obj, theta)
            % Calculate transformation matrix for ply angle theta
            c = cos(theta);
            s = sin(theta);
            T = [c^2, s^2, 2*c*s;
                 s^2, c^2, -2*c*s;
                 -c*s, c*s, c^2-s^2];
        end
        
        function [Me, Ke] = getElementMatrices(obj, xe, ye, ze)
            % Calculate element matrices for composite shell
            % This method is called by ModalAnalysis
            
            % Get element size
            dx = max(xe) - min(xe);
            dy = max(ye) - min(ye);
            
            % Initialize matrices for 6 DOF per node (24x24)
            Me = zeros(24,24);
            Ke = zeros(24,24);
            
            % Gauss quadrature points and weights (2x2 integration)
            xi = [-1/sqrt(3), 1/sqrt(3)];
            w = [1, 1];
            
            % Average density and thickness for mass matrix
            rho_avg = mean(cellfun(@(x) x.rho, obj.Layers));
            
            % Integration loop
            for i = 1:length(xi)
                for j = 1:length(xi)
                    % Current Gauss point
                    xi_i = xi(i);
                    eta_j = xi(j);
                    
                    % Shape functions and derivatives
                    [N, dN] = obj.getShapeFunctions(xi_i, eta_j);
                    
                    % Jacobian matrix
                    J = zeros(2,2);
                    J(1,1) = dx/2;  % dx/dxi
                    J(1,2) = 0;     % dy/dxi
                    J(2,1) = 0;     % dx/deta
                    J(2,2) = dy/2;  % dy/deta
                    detJ = det(J);
                    invJ = inv(J);
                    
                    % Strain-displacement matrix
                    B = zeros(6,24);
                    for k = 1:4
                        % Node DOF indices
                        idx = (k-1)*6 + (1:6);
                        
                        % Transform derivatives to global coordinates
                        dNdx = invJ(1,1)*dN(1,k) + invJ(1,2)*dN(2,k);
                        dNdy = invJ(2,1)*dN(1,k) + invJ(2,2)*dN(2,k);
                        
                        % Fill strain-displacement matrix
                        B(1,idx(1)) = dNdx;        % du/dx
                        B(2,idx(2)) = dNdy;        % dv/dy
                        B(3,idx(3)) = N(k);        % w
                        B(4,idx(4)) = dNdx;        % dθx/dx
                        B(5,idx(5)) = dNdy;        % dθy/dy
                        B(6,idx(6)) = dNdx + dNdy; % dθz/dx + dθz/dy
                    end
                    
                    % Calculate element matrices
                    % Mass matrix
                    for k = 1:4
                        for l = 1:4
                            idx_k = (k-1)*6 + (1:6);
                            idx_l = (l-1)*6 + (1:6);
                            Me(idx_k,idx_l) = Me(idx_k,idx_l) + ...
                                rho_avg*obj.TotalThickness*N(k)*N(l)*detJ*w(i)*w(j);
                        end
                    end
                    
                    % Stiffness matrix using ABD matrix
                    Ke = Ke + B'*obj.ABD*B*detJ*w(i)*w(j);
                end
            end
            
            % Add small regularization for numerical stability
            eps = 1e-6*trace(Ke)/size(Ke,1);
            Ke = Ke + eps*eye(size(Ke));
            Me = Me + eps*eye(size(Me));
        end
        
        function [N, dN] = getShapeFunctions(obj, xi, eta)
            % Calculate shape functions and derivatives for 4-node element
            % xi, eta: Natural coordinates (-1 to 1)
            
            % Shape functions
            N = zeros(1,4);
            N(1) = 0.25*(1-xi)*(1-eta);  % Node 1
            N(2) = 0.25*(1+xi)*(1-eta);  % Node 2
            N(3) = 0.25*(1+xi)*(1+eta);  % Node 3
            N(4) = 0.25*(1-xi)*(1+eta);  % Node 4
            
            % Derivatives with respect to xi and eta
            dN = zeros(2,4);
            % d/dxi
            dN(1,1) = -0.25*(1-eta);
            dN(1,2) = 0.25*(1-eta);
            dN(1,3) = 0.25*(1+eta);
            dN(1,4) = -0.25*(1+eta);
            % d/deta
            dN(2,1) = -0.25*(1-xi);
            dN(2,2) = -0.25*(1+xi);
            dN(2,3) = 0.25*(1+xi);
            dN(2,4) = 0.25*(1-xi);
        end
    end
end
