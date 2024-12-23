# Algorithm Documentation

## Modal Analysis Algorithms

### 1. Subspace Iteration
```matlab
classdef SubspaceIterator < handle
    % Implementation of subspace iteration for modal analysis
    
    properties
        K           % Stiffness matrix
        M           % Mass matrix
        num_modes   % Number of modes to compute
        tolerance   % Convergence tolerance
        max_iter    % Maximum iterations
        block_size  % Size of subspace
    end
    
    methods
        function obj = SubspaceIterator(K, M, num_modes)
            obj.K = K;
            obj.M = M;
            obj.num_modes = num_modes;
            obj.block_size = min(2*num_modes, size(K,1));
            obj.tolerance = 1e-6;
            obj.max_iter = 100;
        end
        
        function [V, D] = iterate(obj)
            n = size(obj.K, 1);
            
            % Initial subspace
            X = rand(n, obj.block_size);
            X = obj.orthogonalize(X);
            
            for iter = 1:obj.max_iter
                % Subspace projection
                KX = obj.K * X;
                MX = obj.M * X;
                
                % Solve reduced eigenproblem
                [Q, Lambda] = eig(X' * KX, X' * MX);
                [lambda, idx] = sort(diag(Lambda));
                Q = Q(:, idx);
                
                % Update subspace
                X_new = X * Q;
                
                % Check convergence
                if obj.checkConvergence(X, X_new)
                    break;
                end
                
                X = X_new;
            end
            
            % Extract results
            V = X(:, 1:obj.num_modes);
            D = diag(lambda(1:obj.num_modes));
        end
        
        function X = orthogonalize(obj, X)
            % M-orthogonalization using modified Gram-Schmidt
            for i = 1:size(X,2)
                % M-normalize current vector
                norm_i = sqrt(X(:,i)' * obj.M * X(:,i));
                X(:,i) = X(:,i) / norm_i;
                
                % M-orthogonalize against remaining vectors
                for j = i+1:size(X,2)
                    proj = X(:,i)' * obj.M * X(:,j);
                    X(:,j) = X(:,j) - proj * X(:,i);
                end
            end
        end
        
        function conv = checkConvergence(obj, X1, X2)
            % Check convergence using angle between subspaces
            angle = subspace(X1, X2);
            conv = angle < obj.tolerance;
        end
    end
end
```

### 2. Block Lanczos
```matlab
classdef BlockLanczos < handle
    % Block Lanczos algorithm for large-scale eigenproblems
    
    properties
        K           % Stiffness matrix
        M           % Mass matrix
        num_modes   % Number of modes to compute
        block_size  % Block size
        max_blocks  % Maximum number of blocks
        tolerance   % Convergence tolerance
    end
    
    methods
        function obj = BlockLanczos(K, M, num_modes)
            obj.K = K;
            obj.M = M;
            obj.num_modes = num_modes;
            obj.block_size = min(num_modes, 10);
            obj.max_blocks = ceil(2*num_modes/obj.block_size);
            obj.tolerance = 1e-6;
        end
        
        function [V, D] = solve(obj)
            n = size(obj.K, 1);
            
            % Initialize
            V0 = zeros(n, obj.block_size);
            V1 = randn(n, obj.block_size);
            V1 = obj.orthogonalize(V1);
            beta = zeros(obj.block_size);
            
            % Storage for Lanczos vectors and coefficients
            V = zeros(n, obj.block_size * obj.max_blocks);
            T = zeros(obj.block_size * obj.max_blocks);
            
            for j = 1:obj.max_blocks
                % Range for current block
                idx = (j-1)*obj.block_size + 1 : j*obj.block_size;
                
                % Store Lanczos vectors
                V(:,idx) = V1;
                
                % Compute K*V1
                KV = obj.K * V1;
                
                % Compute alpha
                alpha = V1' * KV;
                T(idx,idx) = alpha;
                
                if j < obj.max_blocks
                    % Update residual
                    R = KV - V1*alpha - V0*beta';
                    
                    % Compute new beta and V0
                    [V2, beta_new] = qr(R, 0);
                    
                    % Update for next iteration
                    V0 = V1;
                    V1 = V2;
                    beta = beta_new;
                    
                    % Store off-diagonal blocks
                    T(idx,idx+obj.block_size) = beta;
                    T(idx+obj.block_size,idx) = beta';
                end
                
                % Check convergence
                if obj.checkConvergence(T(1:idx,1:idx))
                    break;
                end
            end
            
            % Solve reduced eigenproblem
            [Q, D] = eig(T(1:idx,1:idx));
            [d, p] = sort(diag(D));
            Q = Q(:,p);
            
            % Extract desired eigenpairs
            V = V(:,1:idx) * Q(:,1:obj.num_modes);
            D = diag(d(1:obj.num_modes));
        end
        
        function X = orthogonalize(obj, X)
            % M-orthogonalization
            [Q, R] = qr(obj.M * X, 0);
            X = Q;
        end
        
        function conv = checkConvergence(obj, T)
            % Check convergence using Ritz values
            [~, D] = eig(T);
            d = sort(diag(D));
            if length(d) > obj.num_modes
                conv = abs(d(obj.num_modes+1) - d(obj.num_modes)) < ...
                       obj.tolerance * abs(d(obj.num_modes));
            else
                conv = false;
            end
        end
    end
end
```

## Nonlinear Solution Algorithms

### 1. Arc-Length Method
```matlab
classdef ArcLengthSolver < handle
    % Arc-length method for nonlinear analysis
    
    properties
        K0          % Initial stiffness matrix
        F           % External force vector
        u           % Displacement vector
        lambda      % Load factor
        delta_l     % Arc-length
        max_iter    % Maximum iterations
        tolerance   % Convergence tolerance
    end
    
    methods
        function obj = ArcLengthSolver(K0, F)
            obj.K0 = K0;
            obj.F = F;
            obj.u = zeros(size(F));
            obj.lambda = 0;
            obj.delta_l = 1.0;
            obj.max_iter = 50;
            obj.tolerance = 1e-6;
        end
        
        function [u, lambda] = solve(obj)
            % Predictor step
            [du_pred, dlambda_pred] = obj.predictor();
            
            u = obj.u + du_pred;
            lambda = obj.lambda + dlambda_pred;
            
            for iter = 1:obj.max_iter
                % Calculate residual
                R = obj.lambda*obj.F - obj.internalForce(u);
                
                % Check convergence
                if norm(R) < obj.tolerance
                    break;
                end
                
                % Corrector step
                [du, dlambda] = obj.corrector(R);
                
                % Update solution
                u = u + du;
                lambda = lambda + dlambda;
            end
            
            % Store results
            obj.u = u;
            obj.lambda = lambda;
        end
        
        function [du, dlambda] = predictor(obj)
            % Tangent predictor
            K = obj.tangentStiffness(obj.u);
            du_bar = K \ obj.F;
            
            % Normalize predictor
            du = obj.delta_l * du_bar / norm(du_bar);
            dlambda = obj.delta_l / norm(du_bar);
        end
        
        function [du, dlambda] = corrector(obj, R)
            % Compute correction using cylindrical arc-length
            K = obj.tangentStiffness(obj.u);
            
            % Solve augmented system
            n = length(R);
            A = [K, -obj.F; obj.u', 0];
            b = [R; 0];
            x = A \ b;
            
            du = x(1:n);
            dlambda = x(n+1);
        end
        
        function K = tangentStiffness(obj, u)
            % Calculate tangent stiffness matrix
            % Implementation depends on specific problem
        end
        
        function F_int = internalForce(obj, u)
            % Calculate internal force vector
            % Implementation depends on specific problem
        end
    end
end
```

### 2. Trust Region Algorithm
```matlab
classdef TrustRegionOptimizer < handle
    % Trust region algorithm for optimization
    
    properties
        objective   % Objective function
        gradient    % Gradient function
        hessian    % Hessian function
        x          % Current point
        delta      % Trust region radius
        eta        % Acceptance parameter
        max_iter   % Maximum iterations
        tolerance  % Convergence tolerance
    end
    
    methods
        function obj = TrustRegionOptimizer(objective, gradient, hessian)
            obj.objective = objective;
            obj.gradient = gradient;
            obj.hessian = hessian;
            obj.delta = 1.0;
            obj.eta = 0.1;
            obj.max_iter = 100;
            obj.tolerance = 1e-6;
        end
        
        function [x, f] = optimize(obj, x0)
            obj.x = x0;
            
            for iter = 1:obj.max_iter
                % Compute model
                f = obj.objective(obj.x);
                g = obj.gradient(obj.x);
                H = obj.hessian(obj.x);
                
                % Check convergence
                if norm(g) < obj.tolerance
                    break;
                end
                
                % Solve trust region subproblem
                p = obj.solveSubproblem(g, H);
                
                % Compute actual and predicted reduction
                actual_red = f - obj.objective(obj.x + p);
                pred_red = -(g'*p + 0.5*p'*H*p);
                
                % Update trust region
                rho = actual_red / pred_red;
                if rho > obj.eta
                    obj.x = obj.x + p;
                    if rho > 0.75
                        obj.delta = min(2*obj.delta, 10);
                    end
                else
                    obj.delta = 0.25 * obj.delta;
                end
            end
            
            x = obj.x;
            f = obj.objective(x);
        end
        
        function p = solveSubproblem(obj, g, H)
            % Solve trust region subproblem using Steihaug-Toint CG
            n = length(g);
            p = zeros(n, 1);
            r = g;
            d = -r;
            
            for j = 1:n
                Hd = H*d;
                dHd = d'*Hd;
                
                if dHd <= 0
                    % Negative curvature
                    p = obj.boundaryIntersection(p, d);
                    break;
                end
                
                alpha = (r'*r)/(dHd);
                p_new = p + alpha*d;
                
                if norm(p_new) >= obj.delta
                    % Trust region boundary reached
                    p = obj.boundaryIntersection(p, d);
                    break;
                end
                
                r_new = r + alpha*Hd;
                beta = (r_new'*r_new)/(r'*r);
                d = -r_new + beta*d;
                
                p = p_new;
                r = r_new;
                
                if norm(r) < obj.tolerance
                    break;
                end
            end
        end
        
        function p = boundaryIntersection(obj, p, d)
            % Find intersection with trust region boundary
            a = d'*d;
            b = 2*p'*d;
            c = p'*p - obj.delta^2;
            
            tau = (-b + sqrt(b^2 - 4*a*c))/(2*a);
            p = p + tau*d;
        end
    end
end
```

## Advanced Material Algorithms

### 1. Return Mapping Algorithm
```matlab
classdef ReturnMapping < handle
    % Return mapping algorithm for elastoplastic materials
    
    properties
        E           % Young's modulus
        nu          % Poisson's ratio
        sigma_y     % Yield stress
        H           % Hardening modulus
        tolerance   % Convergence tolerance
        max_iter    % Maximum iterations
    end
    
    methods
        function obj = ReturnMapping(props)
            obj.E = props.E;
            obj.nu = props.nu;
            obj.sigma_y = props.sigma_y;
            obj.H = props.H;
            obj.tolerance = 1e-6;
            obj.max_iter = 50;
        end
        
        function [sigma, ep, alpha] = compute(obj, strain, ep_n, alpha_n)
            % Elastic predictor
            C = obj.elasticityTensor();
            sigma_trial = C * (strain - ep_n);
            s_trial = obj.deviator(sigma_trial);
            q_trial = sqrt(1.5 * sum(s_trial.^2));
            
            % Check yield condition
            f_trial = q_trial - (obj.sigma_y + obj.H*alpha_n);
            
            if f_trial <= 0
                % Elastic step
                sigma = sigma_trial;
                ep = ep_n;
                alpha = alpha_n;
            else
                % Plastic step
                [sigma, ep, alpha] = obj.returnMap(sigma_trial, ep_n, alpha_n);
            end
        end
        
        function [sigma, ep, alpha] = returnMap(obj, sigma_trial, ep_n, alpha_n)
            % Initialize
            dgamma = 0;
            s = obj.deviator(sigma_trial);
            n = sqrt(1.5) * s / norm(s);
            
            for iter = 1:obj.max_iter
                % Update stress and internal variables
                sigma = sigma_trial - 2*obj.G()*dgamma*n;
                alpha = alpha_n + sqrt(2/3)*dgamma;
                
                % Check yield condition
                f = sqrt(1.5*sum(obj.deviator(sigma).^2)) - ...
                    (obj.sigma_y + obj.H*alpha);
                
                if abs(f) < obj.tolerance
                    break;
                end
                
                % Update plastic multiplier
                denom = 2*obj.G() + 2/3*obj.H;
                dgamma = dgamma + f/denom;
            end
            
            % Update plastic strain
            ep = ep_n + dgamma*n;
        end
        
        function C = elasticityTensor(obj)
            % Fourth-order elasticity tensor
            lambda = obj.E*obj.nu/((1+obj.nu)*(1-2*obj.nu));
            mu = obj.E/(2*(1+obj.nu));
            
            C = zeros(6,6);
            C(1:3,1:3) = lambda;
            C = C + 2*mu*eye(6);
        end
        
        function s = deviator(obj, sigma)
            % Compute stress deviator
            p = sum(sigma(1:3))/3;
            s = sigma;
            s(1:3) = s(1:3) - p;
        end
        
        function G = G(obj)
            % Shear modulus
            G = obj.E/(2*(1+obj.nu));
        end
    end
end
```

## Advanced Numerical Methods

### 1. Adaptive Time Integration
```matlab
classdef AdaptiveTimeIntegrator < handle
    % Adaptive time integration with error control
    
    properties
        M           % Mass matrix
        C           % Damping matrix
        K           % Stiffness matrix
        F           % External force
        tolerance   % Error tolerance
        dt_min      % Minimum time step
        dt_max      % Maximum time step
    end
    
    methods
        function obj = AdaptiveTimeIntegrator(M, C, K, F)
            obj.M = M;
            obj.C = C;
            obj.K = K;
            obj.F = F;
            obj.tolerance = 1e-6;
            obj.dt_min = 1e-6;
            obj.dt_max = 1e-2;
        end
        
        function [t, u, v, a] = solve(obj, tspan)
            % Initialize
            t = tspan(1);
            dt = obj.dt_max;
            u = zeros(size(obj.M,1), 1);
            v = zeros(size(u));
            a = obj.M \ (obj.F - obj.C*v - obj.K*u);
            
            % Storage
            t_hist = t;
            u_hist = u;
            v_hist = v;
            a_hist = a;
            
            while t < tspan(2)
                % Predict with current step
                [u1, v1, a1] = obj.step(u, v, a, dt);
                
                % Predict with half steps
                dt_half = dt/2;
                [u_half, v_half, a_half] = obj.step(u, v, a, dt_half);
                [u2, v2, a2] = obj.step(u_half, v_half, a_half, dt_half);
                
                % Estimate error
                error = norm(u2 - u1) / norm(u2);
                
                % Adjust time step
                if error > obj.tolerance
                    % Reduce time step
                    dt = max(dt/2, obj.dt_min);
                    continue;
                else
                    % Accept step and possibly increase dt
                    t = t + dt;
                    u = u2;
                    v = v2;
                    a = a2;
                    
                    % Store results
                    t_hist = [t_hist; t];
                    u_hist = [u_hist, u];
                    v_hist = [v_hist, v];
                    a_hist = [a_hist, a];
                    
                    if error < obj.tolerance/10
                        dt = min(2*dt, obj.dt_max);
                    end
                end
            end
            
            % Output history
            t = t_hist;
            u = u_hist;
            v = v_hist;
            a = a_hist;
        end
        
        function [u_new, v_new, a_new] = step(obj, u, v, a, dt)
            % Newmark-beta method with parameters for unconditional stability
            beta = 0.25;
            gamma = 0.5;
            
            % Predict
            u_pred = u + dt*v + dt^2/2*((1-2*beta)*a);
            v_pred = v + dt*(1-gamma)*a;
            
            % Effective stiffness
            K_eff = obj.K + obj.M/(beta*dt^2) + obj.C*gamma/(beta*dt);
            
            % Effective force
            F_eff = obj.F - obj.K*u_pred - obj.C*v_pred;
            
            % Solve for acceleration
            a_new = K_eff \ F_eff;
            
            % Correct
            u_new = u_pred + beta*dt^2*a_new;
            v_new = v_pred + gamma*dt*a_new;
        end
    end
end
```

### 2. Multi-Scale Homogenization
```matlab
classdef MicroScaleHomogenizer < handle
    % Microscale homogenization for composite materials
    
    properties
        rve_size        % RVE dimensions
        mesh            % Microscale mesh
        fiber_props     % Fiber properties
        matrix_props    % Matrix properties
        damage_model    % Microscale damage model
    end
    
    methods
        function obj = MicroScaleHomogenizer(rve_params)
            obj.rve_size = rve_params.size;
            obj.fiber_props = rve_params.fiber;
            obj.matrix_props = rve_params.matrix;
            obj.createMesh(rve_params.mesh_size);
        end
        
        function createMesh(obj, mesh_size)
            % Generate RVE mesh with periodic boundaries
            obj.mesh = obj.generatePeriodicMesh(mesh_size);
            obj.identifyPeriodic();
        end
        
        function [C_eff, stress] = homogenize(obj, strain)
            % Compute effective properties through homogenization
            
            % Apply periodic boundary conditions
            [K, dofs] = obj.assembleStiffness();
            u_bc = obj.applyPeriodicBC(strain);
            
            % Solve microscale problem
            u = obj.solveMicroscale(K, dofs, u_bc);
            
            % Compute microscale stresses
            stress = obj.computeMicroStress(u);
            
            % Homogenize properties
            C_eff = obj.computeEffectiveProperties(stress, strain);
        end
        
        function updateDamage(obj, stress)
            % Update microscale damage based on stress state
            if ~isempty(obj.damage_model)
                damage = obj.damage_model.evaluate(stress);
                obj.updateMaterialProperties(damage);
            end
        end
        
        function [K, dofs] = assembleStiffness(obj)
            % Assemble microscale stiffness with damage
            K = sparse(obj.mesh.ndof, obj.mesh.ndof);
            
            for e = 1:obj.mesh.nelem
                % Element properties with damage
                props = obj.getElementProperties(e);
                
                % Element stiffness
                Ke = obj.elementStiffness(props);
                
                % Assembly
                dofs = obj.mesh.elementDofs(e);
                K(dofs, dofs) = K(dofs, dofs) + Ke;
            end
            
            dofs = 1:obj.mesh.ndof;
        end
        
        function u = solveMicroscale(obj, K, dofs, u_bc)
            % Solve microscale boundary value problem
            
            % Apply boundary conditions
            [K_mod, f_mod] = obj.modifySystem(K, dofs, u_bc);
            
            % Solve system
            u = K_mod \ f_mod;
            
            % Recover full solution
            u = obj.recoverSolution(u, u_bc);
        end
        
        function C_eff = computeEffectiveProperties(obj, stress, strain)
            % Compute effective properties through volume averaging
            
            % Volume average of stress and strain
            stress_avg = obj.volumeAverage(stress);
            strain_avg = obj.volumeAverage(strain);
            
            % Effective stiffness
            C_eff = stress_avg / strain_avg;
        end
        
        function stress = computeMicroStress(obj, u)
            % Compute microscale stress field
            stress = zeros(obj.mesh.nelem, 6);
            
            for e = 1:obj.mesh.nelem
                % Element strain
                strain_e = obj.elementStrain(u, e);
                
                % Element properties
                props = obj.getElementProperties(e);
                
                % Element stress
                stress(e,:) = props.C * strain_e;
            end
        end
        
        function avg = volumeAverage(obj, field)
            % Compute volume average of field
            V_total = 0;
            avg = zeros(size(field,2), 1);
            
            for e = 1:obj.mesh.nelem
                V_e = obj.mesh.elementVolume(e);
                avg = avg + V_e * field(e,:)';
                V_total = V_total + V_e;
            end
            
            avg = avg / V_total;
        end
    end
end
```

### 3. Contact Detection Algorithm
```matlab
classdef ContactDetector < handle
    % Contact detection and force computation
    
    properties
        master      % Master surface
        slave       % Slave surface
        gap_tol     % Contact tolerance
        friction    % Friction coefficient
        penalty     % Normal penalty parameter
        tang_penalty % Tangential penalty
    end
    
    methods
        function obj = ContactDetector(master, slave)
            obj.master = master;
            obj.slave = slave;
            obj.gap_tol = 1e-6;
            obj.friction = 0;
            obj.penalty = 1e6;
            obj.tang_penalty = 1e5;
        end
        
        function [forces, status] = detectContact(obj, x_master, x_slave)
            % Detect contact and compute forces
            
            forces = zeros(size(x_slave));
            status = zeros(size(x_slave,1), 1);
            
            % Build spatial search structure
            tree = obj.buildSearchTree(x_master);
            
            % Check each slave node
            for i = 1:size(x_slave,1)
                % Find closest master segment
                [seg_id, proj] = obj.findClosestSegment(tree, x_slave(i,:));
                
                % Compute gap
                gap = norm(x_slave(i,:) - proj);
                
                if gap < obj.gap_tol
                    % Contact detected
                    status(i) = 1;
                    
                    % Normal force
                    n = (x_slave(i,:) - proj) / gap;
                    f_n = obj.penalty * gap * n;
                    
                    % Tangential force (friction)
                    if obj.friction > 0
                        v_t = obj.relativeTangentialVelocity(i, seg_id);
                        f_t = obj.computeFrictionForce(v_t, f_n);
                        forces(i,:) = f_n + f_t;
                    else
                        forces(i,:) = f_n;
                    end
                end
            end
        end
        
        function tree = buildSearchTree(obj, x)
            % Build k-d tree for spatial search
            tree = KDTree(x);
        end
        
        function [seg_id, proj] = findClosestSegment(obj, tree, x)
            % Find closest master segment and projection
            
            % Find candidate segments
            [idx, dist] = tree.findKNearest(x, 10);
            
            % Check each candidate
            min_dist = inf;
            seg_id = 0;
            proj = zeros(1,3);
            
            for i = 1:length(idx)
                % Get segment nodes
                seg = obj.master.getSegment(idx(i));
                
                % Project point
                [p, d] = obj.projectPointToSegment(x, seg);
                
                if d < min_dist
                    min_dist = d;
                    seg_id = idx(i);
                    proj = p;
                end
            end
        end
        
        function [p, d] = projectPointToSegment(obj, x, seg)
            % Project point onto segment
            
            % Segment vectors
            v1 = seg(2,:) - seg(1,:);
            v2 = x - seg(1,:);
            
            % Projection parameter
            t = dot(v2,v1) / dot(v1,v1);
            t = max(0, min(1, t));
            
            % Projection point
            p = seg(1,:) + t*v1;
            
            % Distance
            d = norm(x - p);
        end
        
        function v_t = relativeTangentialVelocity(obj, slave_id, master_id)
            % Compute relative tangential velocity
            v_s = obj.slave.getVelocity(slave_id);
            v_m = obj.master.getVelocity(master_id);
            n = obj.master.getNormal(master_id);
            
            % Remove normal component
            v_rel = v_s - v_m;
            v_t = v_rel - dot(v_rel,n)*n;
        end
        
        function f_t = computeFrictionForce(obj, v_t, f_n)
            % Compute friction force
            v_mag = norm(v_t);
            
            if v_mag > 0
                % Sliding friction
                f_t = -obj.friction * norm(f_n) * v_t/v_mag;
            else
                % Static friction (regularized)
                f_t = -obj.tang_penalty * v_t;
            end
        end
    end
end
