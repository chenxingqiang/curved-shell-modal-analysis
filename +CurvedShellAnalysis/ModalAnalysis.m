classdef ModalAnalysis < handle
    % ModalAnalysis Class for modal analysis of curved shells
    
    properties
        Surface         % Surface object
        NumModes        % Number of modes to compute
        frequencies     % Natural frequencies
        modes          % Mode shapes
        mesh           % Mesh data
    end
    
    methods
        function obj = ModalAnalysis(surface, num_modes)
            % Constructor for ModalAnalysis
            obj.Surface = surface;
            obj.NumModes = num_modes;
        end
        
        function analyze(obj)
            % Perform modal analysis
            fprintf('开始模态分析...\n');
            
            % Get mesh from surface
            obj.mesh = obj.Surface.getMesh();
            
            % Get material properties
            E = obj.Surface.E;
            nu = obj.Surface.nu;
            rho = obj.Surface.rho;
            t = obj.Surface.t;
            
            % Get mesh dimensions
            [ny, nx] = size(obj.mesh.X);
            ndof = 3;  % 3 DOF per node (u, v, w)
            
            % Initialize matrices
            n = nx * ny * ndof;
            M = sparse(n, n);
            K = sparse(n, n);
            
            % Assembly
            fprintf('组装质量和刚度矩阵...\n');
            for i = 1:nx-1
                for j = 1:ny-1
                    % Element nodes
                    nodes = [j+(i-1)*ny, j+1+(i-1)*ny, j+1+i*ny, j+i*ny];
                    
                    % Element coordinates
                    xe = obj.mesh.X(nodes([1,2,3,4]));
                    ye = obj.mesh.Y(nodes([1,2,3,4]));
                    ze = obj.mesh.Z(nodes([1,2,3,4]));
                    
                    % Element matrices
                    [Me, Ke] = obj.elementMatrices(xe, ye, ze, E, nu, rho, t);
                    
                    % Global DOFs
                    dofs = [];
                    for node = nodes
                        dofs = [dofs, (node-1)*ndof+1:node*ndof];
                    end
                    
                    % Assembly
                    M(dofs,dofs) = M(dofs,dofs) + Me;
                    K(dofs,dofs) = K(dofs,dofs) + Ke;
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
            
            % Solve eigenvalue problem
            fprintf('求解特征值问题...\n');
            
            % Add small stiffness to prevent singularity
            K = K + sparse(1:n, 1:n, 1e-6*max(diag(K)), n, n);
            
            % Solve reduced system
            Kr = K(free_dofs,free_dofs);
            Mr = M(free_dofs,free_dofs);
            
            % Use shift-invert mode for better numerical stability
            sigma = 1e-6;  % Small shift
            [V, D] = eigs(Kr, Mr, obj.NumModes, sigma);
            
            % Store results
            obj.frequencies = sqrt(diag(D))/(2*pi);
            obj.modes = zeros(n, obj.NumModes);
            obj.modes(free_dofs,:) = V;
        end
        
        function response = calculateTransientResponse(obj, F, t)
            % Calculate transient response to force history F at times t
            % F: Force history vector
            % t: Time vector
            
            % Default damping ratio
            zeta = 0.02;
            
            % Initialize response
            response = zeros(size(t));
            
            % Modal participation factors
            for i = 1:obj.NumModes
                % Natural frequency in rad/s
                wn = 2*pi*obj.frequencies(i);
                
                % Damped natural frequency
                wd = wn*sqrt(1-zeta^2);
                
                % Modal response using convolution integral
                h = exp(-zeta*wn*t).*sin(wd*t)./(wd);  % Impulse response
                y = conv(F, h);  % Convolution
                y = y(1:length(t))*mean(diff(t));  % Truncate and scale
                
                % Add modal contribution
                response = response + y;
            end
        end
        
        function energy = calculateModalEnergy(obj, mode_num)
            % Calculate modal strain energy distribution for a given mode
            % mode_num: Mode number to analyze
            
            % Get mesh dimensions
            [ny, nx] = size(obj.mesh.X);
            ndof = 3;  % 3 DOF per node (u, v, w)
            
            % Get material properties
            E = obj.Surface.E;
            nu = obj.Surface.nu;
            t = obj.Surface.t;
            
            % Initialize energy field
            energy = zeros(ny, nx);
            
            % Calculate strain energy at each element
            for i = 1:nx-1
                for j = 1:ny-1
                    % Element nodes
                    nodes = [j+(i-1)*ny, j+1+(i-1)*ny, j+1+i*ny, j+i*ny];
                    
                    % Element coordinates
                    xe = obj.mesh.X(nodes([1,2,3,4]));
                    ye = obj.mesh.Y(nodes([1,2,3,4]));
                    ze = obj.mesh.Z(nodes([1,2,3,4]));
                    
                    % Element displacements
                    dofs = [];
                    for node = nodes
                        dofs = [dofs, (node-1)*ndof+1:node*ndof];
                    end
                    u_e = obj.modes(dofs,mode_num);
                    
                    % Calculate element strain energy
                    [~, Ke] = obj.elementMatrices(xe, ye, ze, E, nu, 1, t);  % Use unit density
                    e = 0.5 * u_e' * Ke * u_e;
                    
                    % Distribute energy to nodes
                    energy(nodes) = energy(nodes) + e/4;  % Average to nodes
                end
            end
        end
        
        function mac = calculateMAC(obj, num_modes)
            % Calculate Modal Assurance Criterion (MAC) matrix
            % num_modes: Number of modes to include in MAC calculation
            
            mac = zeros(num_modes);
            for i = 1:num_modes
                for j = 1:num_modes
                    % Get mode shapes
                    phi_i = obj.modes(:,i);
                    phi_j = obj.modes(:,j);
                    
                    % Calculate MAC value
                    mac(i,j) = abs(phi_i' * phi_j)^2 / ...
                              ((phi_i' * phi_i) * (phi_j' * phi_j));
                end
            end
        end
        
        function plotMode(obj, mode_num)
            % Plot mode shape
            if nargin < 2
                mode_num = 1;
            end
            
            % Get mesh dimensions
            [ny, nx] = size(obj.mesh.X);
            ndof = 3;
            
            % Extract mode shape
            mode = reshape(obj.modes(:,mode_num), ndof, []);
            w = reshape(mode(3,:), ny, nx);
            
            % Plot
            surf(obj.mesh.X, obj.mesh.Y, obj.mesh.Z + w*0.1*max(abs(obj.mesh.Z(:))));
            shading interp;
            colormap('jet');
            axis equal;
            view(45, 30);
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
        end
        
        function [Me, Ke] = elementMatrices(obj, xe, ye, ze, E, nu, rho, t)
            % Calculate element matrices using simplified shell theory
            
            % Element size
            dx = diff(xe([1,2]));
            dy = diff(ye([1,4]));
            
            % Material matrix (plane stress)
            c = E/(1-nu^2);
            D = zeros(3,3);
            D(1:2,1:2) = c * [1, nu; nu, 1];
            D(3,3) = c * (1-nu)/2;  % Shear modulus
            
            % Shape functions at Gauss points
            xi = [-1/sqrt(3), 1/sqrt(3)];
            eta = [-1/sqrt(3), 1/sqrt(3)];
            w = [1, 1];  % Gauss weights
            
            % Initialize element matrices
            Me = zeros(12,12);
            Ke = zeros(12,12);
            
            % Node indices for shape functions
            xi_nodes = [-1, 1, 1, -1];
            eta_nodes = [-1, -1, 1, 1];
            
            % Gauss integration
            for i = 1:2
                for j = 1:2
                    % Shape function derivatives
                    dN = zeros(2,4);
                    dN(1,:) = [-(1-eta(j)), (1-eta(j)), (1+eta(j)), -(1+eta(j))]/4;
                    dN(2,:) = [-(1-xi(i)), -(1+xi(i)), (1+xi(i)), (1-xi(i))]/4;
                    
                    % Jacobian
                    J = zeros(2,2);
                    for k = 1:4
                        J = J + [dN(1,k)*xe(k), dN(1,k)*ye(k);
                                dN(2,k)*xe(k), dN(2,k)*ye(k)];
                    end
                    detJ = det(J);
                    invJ = inv(J);
                    
                    % B matrix
                    B = zeros(3,12);
                    for k = 1:4
                        % Transform derivatives to global coordinates
                        dNdx = invJ(1,1)*dN(1,k) + invJ(1,2)*dN(2,k);
                        dNdy = invJ(2,1)*dN(1,k) + invJ(2,2)*dN(2,k);
                        
                        % Fill B matrix
                        idx = (k-1)*3 + (1:3);
                        B(:,idx) = [dNdx,    0,  0;
                                  0,    dNdy,  0;
                                  dNdy,  dNdx,  0];
                    end
                    
                    % Shape functions for mass matrix
                    N = zeros(3,12);
                    for k = 1:4
                        N_k = (1 + xi(i)*xi_nodes(k))*(1 + eta(j)*eta_nodes(k))/4;
                        idx = (k-1)*3 + (1:3);
                        N(:,idx) = eye(3)*N_k;
                    end
                    
                    % Mass matrix contribution
                    Me = Me + rho*t*N'*N*detJ*w(i)*w(j);
                    
                    % Stiffness matrix contribution
                    Ke = Ke + t*B'*D*B*detJ*w(i)*w(j);
                end
            end
            
            % Add bending stiffness
            D_bend = E*t^3/(12*(1-nu^2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
            for k = 1:4
                idx = (k-1)*3 + 3;  % w-DOF
                Ke(idx,idx) = Ke(idx,idx) + D_bend(1,1);
            end
            
            % Add small regularization term to improve numerical stability
            eps = 1e-10;
            Ke = Ke + eps*eye(size(Ke));
        end
    end
end
