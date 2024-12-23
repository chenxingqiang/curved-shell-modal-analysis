classdef ModalAnalysis < handle
    % ModalAnalysis Class for modal analysis of curved shells
    %   This class performs modal analysis on curved shells, including
    %   composite shells with multiple layers
    
    properties
        Surface         % Surface object
        NumModes        % Number of modes to compute
        frequencies     % Natural frequencies
        modes          % Mode shapes
        mesh           % Mesh data
        damping_ratio  % Modal damping ratio
    end
    
    methods
        function obj = ModalAnalysis(surface, num_modes)
            % Constructor for ModalAnalysis
            obj.Surface = surface;
            obj.NumModes = num_modes;
            if isprop(surface, 'zeta')
                obj.damping_ratio = surface.zeta;
            else
                obj.damping_ratio = 0.02;  % Default damping ratio
            end
        end
        
        function analyze(obj)
            % Perform modal analysis
            fprintf('开始模态分析...\n');
            
            % Get mesh from surface
            obj.mesh = obj.Surface.getMesh();
            
            % Get mesh dimensions
            [ny, nx] = size(obj.mesh.X);
            ndof = 6;  % 6 DOF per node for shell elements
            
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
                    if isa(obj.Surface, 'CurvedShellAnalysis.CompositeShell')
                        [Me, Ke] = obj.Surface.getElementMatrices(xe, ye, ze);
                    else
                        % For non-composite shells, use standard element matrices
                        [Me, Ke] = obj.elementMatrices(xe, ye, ze, ...
                            obj.Surface.E, obj.Surface.nu, obj.Surface.rho, obj.Surface.t);
                    end
                    
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
            
            % Fix all edges
            fixed_nodes = unique([
                1:ny,                    % Bottom edge
                1:ny:ny*(nx-1)+1,       % Left edge
                ny:ny:ny*nx,            % Right edge
                ny*(nx-1)+1:ny*nx       % Top edge
            ]);
            
            % Initialize fixed DOFs array
            fixed_dofs = [];
            
            % For each fixed node, constrain all DOFs
            for node = fixed_nodes
                fixed_dofs = [fixed_dofs, (node-1)*ndof + (1:ndof)];
            end
            
            % Get free DOFs
            free_dofs = setdiff(1:n, fixed_dofs);
            
            % Extract free DOF matrices
            Kr = K(free_dofs,free_dofs);
            Mr = M(free_dofs,free_dofs);
            
            % Add small regularization
            eps = 1e-10;
            Kr = Kr + eps*speye(size(Kr));
            Mr = Mr + eps*speye(size(Mr));
            
            % Scale matrices for better conditioning
            scale = max(abs([Mr(:); Kr(:)]));
            if scale > 0
                Mr = Mr/scale;
                Kr = Kr/scale;
            end
            
            % Solve generalized eigenvalue problem using eigs
            opts.disp = 0;
            opts.isreal = true;
            opts.tol = 1e-6;
            opts.maxit = 1000;
            
            try
                % Try shift-invert mode first
                sigma = 1e-6;  % Small shift
                [V, D] = eigs(Kr, Mr, min(obj.NumModes*2, size(Kr,1)), sigma, opts);
            catch
                warning('Shift-invert mode failed, trying regular mode...');
                try
                    [V, D] = eigs(Kr, Mr, min(obj.NumModes*2, size(Kr,1)), 'sm', opts);
                catch
                    warning('eigs failed, trying full eigenvalue decomposition...');
                    [V, D] = eig(full(Kr), full(Mr));
                end
            end
            
            % Extract eigenvalues and sort
            [lambda, idx] = sort(real(diag(D)), 'ascend');
            
            % Remove negative or complex eigenvalues
            valid = find(lambda > 0 & abs(imag(lambda)./real(lambda)) < 1e-6);
            if isempty(valid)
                error('No valid eigenvalues found. Check model parameters.');
            end
            
            lambda = lambda(valid);
            V = V(:,idx(valid));
            
            % Limit to requested number of modes
            n_modes = min(obj.NumModes, length(lambda));
            lambda = lambda(1:n_modes);
            V = V(:,1:n_modes);
            
            % Update number of modes
            obj.NumModes = n_modes;
            
            % Convert eigenvalues to frequencies (Hz)
            obj.frequencies = sqrt(abs(lambda))/(2*pi)*sqrt(scale);
            
            % Store mode shapes
            obj.modes = zeros(n, n_modes);
            obj.modes(free_dofs,:) = V;
            
            % Normalize mode shapes
            for i = 1:n_modes
                % Mass normalization
                norm_factor = sqrt(obj.modes(:,i)'*M*obj.modes(:,i));
                if norm_factor > 0
                    obj.modes(:,i) = obj.modes(:,i)/norm_factor;
                end
            end
            
            fprintf('找到 %d 个模态，频率范围: %.1f - %.1f Hz\n', ...
                   n_modes, full(min(obj.frequencies)), full(max(obj.frequencies)));
        end
        
        function response = calculateTransientResponse(obj, F, t)
            % Calculate transient response to force history F at times t
            % F: Force history vector
            % t: Time vector
            
            % Modal transformation
            modal_force = obj.modes' * F;
            omega = 2*pi*obj.frequencies;
            
            % Initialize response
            response = zeros(size(t));
            
            % Calculate response for each mode
            for i = 1:obj.NumModes
                % Modal parameters
                wn = omega(i);
                zeta = obj.damping_ratio;
                wd = wn*sqrt(1-zeta^2);  % Damped natural frequency
                
                % Modal response
                h = exp(-zeta*wn*t).*(cos(wd*t) + zeta*wn/wd*sin(wd*t));
                response = response + modal_force(i)*h;
            end
        end
        
        function plotMode(obj, mode_num)
            % Plot mode shape
            if mode_num > obj.NumModes
                error('Mode number exceeds available modes');
            end
            
            % Get mesh dimensions
            [ny, nx] = size(obj.mesh.X);
            n_nodes = nx * ny;
            ndof = size(obj.modes,1)/n_nodes;
            
            % Extract displacement components
            mode = obj.modes(:,mode_num);
            
            % Reshape mode to get out-of-plane displacement
            w = zeros(ny, nx);
            for i = 1:nx
                for j = 1:ny
                    node = j + (i-1)*ny;
                    w(j,i) = mode((node-1)*ndof + 3);  % w displacement is 3rd DOF
                end
            end
            
            % Create deformed surface
            Z = obj.mesh.Z + real(w);
            
            % Plot
            surf(obj.mesh.X, obj.mesh.Y, Z);
            axis equal;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            colormap('jet');
            shading interp;
            
            % Add frequency to title
            title(sprintf('Mode %d: %.1f Hz', mode_num, full(obj.frequencies(mode_num))));
        end
        
        function energy = calculateModalEnergy(obj, mode_num)
            % Calculate modal strain energy distribution for a given mode
            % mode_num: Mode number to analyze
            
            % Get mesh dimensions
            [ny, nx] = size(obj.mesh.X);
            ndof = size(obj.modes,1)/(nx*ny);
            
            % Initialize energy field
            energy = zeros(ny, nx);
            
            % Get mode shape
            mode = obj.modes(:,mode_num);
            
            % Calculate energy for each element
            for i = 1:nx-1
                for j = 1:ny-1
                    % Element nodes
                    nodes = [j+(i-1)*ny, j+1+(i-1)*ny, j+1+i*ny, j+i*ny];
                    
                    % Element coordinates
                    xe = obj.mesh.X(nodes([1,2,3,4]));
                    ye = obj.mesh.Y(nodes([1,2,3,4]));
                    ze = obj.mesh.Z(nodes([1,2,3,4]));
                    
                    % Element DOFs
                    dofs = [];
                    for node = nodes
                        dofs = [dofs, (node-1)*ndof+1:node*ndof];
                    end
                    
                    % Element mode shape
                    q = mode(dofs);
                    
                    % Element matrices
                    if isa(obj.Surface, 'CurvedShellAnalysis.CompositeShell')
                        [~, Ke] = obj.Surface.getElementMatrices(xe, ye, ze);
                    else
                        [~, Ke] = obj.elementMatrices(xe, ye, ze, ...
                            obj.Surface.E, obj.Surface.nu, obj.Surface.rho, obj.Surface.t);
                    end
                    
                    % Element strain energy
                    e = 0.5 * q' * Ke * q;
                    
                    % Distribute energy to nodes
                    energy(nodes) = energy(nodes) + e/4;
                end
            end
        end
    end
end
