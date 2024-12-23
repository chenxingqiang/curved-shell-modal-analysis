classdef ShapeOptimization < handle
    % SHAPEOPTIMIZATION Class for shape and thickness optimization
    
    properties
        Surface  % Surface object
        Analysis  % Analysis object (Modal, Nonlinear, etc.)
        ObjectiveType  % Type of objective ('frequency', 'weight', 'stress')
        Constraints  % Optimization constraints
        DesignVariables  % Design variables and their bounds
        Algorithm  % Optimization algorithm ('genetic', 'gradient', 'pattern')
        MaxIterations = 100  % Maximum optimization iterations
        Tolerance = 1e-6  % Convergence tolerance
        OutputFolder = 'OptimizationResults'  % Output folder
    end
    
    methods
        function obj = ShapeOptimization(surface, analysis)
            % Constructor
            obj.Surface = surface;
            obj.Analysis = analysis;
            
            % Create output folder
            if ~exist(obj.OutputFolder, 'dir')
                mkdir(obj.OutputFolder);
            end
        end
        
        function setObjective(obj, type, params)
            % Set optimization objective
            obj.ObjectiveType = type;
            
            switch type
                case 'frequency'
                    % Maximize specific frequencies or frequency gaps
                    obj.ObjectiveType = struct('type', type, ...
                                             'mode_indices', params.modes, ...
                                             'weights', params.weights);
                case 'weight'
                    % Minimize weight with stress/displacement constraints
                    obj.ObjectiveType = struct('type', type, ...
                                             'stress_limit', params.stress_limit, ...
                                             'disp_limit', params.disp_limit);
                case 'stress'
                    % Minimize maximum stress
                    obj.ObjectiveType = struct('type', type, ...
                                             'weight_limit', params.weight_limit);
            end
        end
        
        function addConstraint(obj, type, params)
            % Add optimization constraint
            if isempty(obj.Constraints)
                obj.Constraints = {};
            end
            
            constraint = struct('type', type, 'params', params);
            obj.Constraints{end+1} = constraint;
        end
        
        function setDesignVariables(obj, variables)
            % Set design variables and their bounds
            obj.DesignVariables = variables;
            
            % Initialize design vector
            obj.DesignVariables.x0 = [];
            obj.DesignVariables.lb = [];
            obj.DesignVariables.ub = [];
            
            % Process each variable
            for i = 1:length(variables)
                var = variables{i};
                switch var.type
                    case 'thickness'
                        obj.DesignVariables.x0(end+1) = obj.Surface.t;
                        obj.DesignVariables.lb(end+1) = var.bounds(1);
                        obj.DesignVariables.ub(end+1) = var.bounds(2);
                    case 'radius'
                        obj.DesignVariables.x0(end+1) = obj.Surface.R;
                        obj.DesignVariables.lb(end+1) = var.bounds(1);
                        obj.DesignVariables.ub(end+1) = var.bounds(2);
                    case 'shape'
                        n_params = length(var.control_points);
                        obj.DesignVariables.x0(end+1:end+n_params) = var.initial_values;
                        obj.DesignVariables.lb(end+1:end+n_params) = var.bounds(1);
                        obj.DesignVariables.ub(end+1:end+n_params) = var.bounds(2);
                end
            end
        end
        
        function [x_opt, f_opt, history] = optimize(obj)
            % Perform optimization
            fprintf('开始优化...\n');
            
            % Set optimization options
            options = optimoptions('fmincon', ...
                'Display', 'iter', ...
                'MaxIterations', obj.MaxIterations, ...
                'OptimalityTolerance', obj.Tolerance, ...
                'UseParallel', true);
            
            % Define objective and constraint functions
            obj_fun = @(x) obj.evaluateObjective(x);
            con_fun = @(x) obj.evaluateConstraints(x);
            
            % Run optimization
            [x_opt, f_opt, ~, output] = fmincon(obj_fun, ...
                obj.DesignVariables.x0, ...
                [], [], [], [], ...
                obj.DesignVariables.lb, ...
                obj.DesignVariables.ub, ...
                con_fun, options);
            
            % Store optimization history
            history = struct('iterations', output.iterations, ...
                           'funcCount', output.funcCount, ...
                           'firstorderopt', output.firstorderopt);
            
            % Update surface with optimal design
            obj.updateDesign(x_opt);
            
            % Plot results
            obj.plotOptimizationResults(history);
        end
        
        function f = evaluateObjective(obj, x)
            % Evaluate objective function
            
            % Update design
            obj.updateDesign(x);
            
            % Calculate objective based on type
            switch obj.ObjectiveType.type
                case 'frequency'
                    % Run modal analysis
                    obj.Analysis.analyze();
                    freqs = obj.Analysis.NaturalFrequencies;
                    
                    % Calculate weighted sum of frequencies
                    f = -sum(obj.ObjectiveType.weights .* freqs(obj.ObjectiveType.mode_indices));
                    
                case 'weight'
                    % Calculate total weight
                    [X, Y, Z] = obj.Surface.generateMesh();
                    area = sum(sum(sqrt(gradient(X).^2 + gradient(Y).^2 + gradient(Z).^2)));
                    f = area * obj.Surface.t * obj.Surface.rho;
                    
                case 'stress'
                    % Run stress analysis
                    if isa(obj.Analysis, 'NonlinearAnalysis')
                        [~, stress, ~] = obj.Analysis.solveNonlinear(obj.Analysis.LoadVector);
                    else
                        stress = obj.Analysis.calculateStress();
                    end
                    f = max(max(stress.vm));
            end
        end
        
        function [c, ceq] = evaluateConstraints(obj, x)
            % Evaluate constraints
            
            % Update design
            obj.updateDesign(x);
            
            % Initialize constraint vectors
            c = [];
            ceq = [];
            
            % Evaluate each constraint
            for i = 1:length(obj.Constraints)
                constraint = obj.Constraints{i};
                
                switch constraint.type
                    case 'stress'
                        % Stress constraint
                        if isa(obj.Analysis, 'NonlinearAnalysis')
                            [~, stress, ~] = obj.Analysis.solveNonlinear(obj.Analysis.LoadVector);
                        else
                            stress = obj.Analysis.calculateStress();
                        end
                        c(end+1) = max(max(stress.vm)) - constraint.params.limit;
                        
                    case 'displacement'
                        % Displacement constraint
                        if isa(obj.Analysis, 'NonlinearAnalysis')
                            [U, ~, ~] = obj.Analysis.solveNonlinear(obj.Analysis.LoadVector);
                        else
                            U = obj.Analysis.solve();
                        end
                        c(end+1) = max(abs(U)) - constraint.params.limit;
                        
                    case 'frequency'
                        % Frequency constraint
                        obj.Analysis.analyze();
                        freqs = obj.Analysis.NaturalFrequencies;
                        c(end+1) = constraint.params.min_freq - min(freqs);
                        
                    case 'volume'
                        % Volume constraint
                        [X, Y, Z] = obj.Surface.generateMesh();
                        vol = sum(sum(sqrt(gradient(X).^2 + gradient(Y).^2 + gradient(Z).^2))) * obj.Surface.t;
                        c(end+1) = vol - constraint.params.limit;
                end
            end
        end
        
        function updateDesign(obj, x)
            % Update design variables in the surface object
            idx = 1;
            
            for i = 1:length(obj.DesignVariables)
                var = obj.DesignVariables{i};
                switch var.type
                    case 'thickness'
                        obj.Surface.t = x(idx);
                        idx = idx + 1;
                    case 'radius'
                        obj.Surface.R = x(idx);
                        idx = idx + 1;
                    case 'shape'
                        n_params = length(var.control_points);
                        obj.Surface.updateShape(var.control_points, x(idx:idx+n_params-1));
                        idx = idx + n_params;
                end
            end
        end
        
        function plotOptimizationResults(obj, history)
            % Plot optimization results
            
            % Convergence history
            figure('Position', [100 100 800 600]);
            subplot(2,2,1);
            plot(1:history.iterations, history.firstorderopt, 'b-o', 'LineWidth', 2);
            grid on;
            title('优化收敛历史');
            xlabel('迭代次数');
            ylabel('一阶最优性');
            
            % Original vs optimized shape
            subplot(2,2,2);
            [X, Y, Z] = obj.Surface.generateMesh();
            surf(X, Y, Z);
            title('优化后形状');
            axis equal;
            
            % Objective function value
            subplot(2,2,3);
            switch obj.ObjectiveType.type
                case 'frequency'
                    obj.Analysis.analyze();
                    freqs = obj.Analysis.NaturalFrequencies;
                    bar(obj.ObjectiveType.mode_indices, freqs(obj.ObjectiveType.mode_indices));
                    title('优化后固有频率');
                    xlabel('模态序号');
                    ylabel('频率 (Hz)');
                    
                case 'weight'
                    text(0.5, 0.5, sprintf('优化后重量: %.2f kg', obj.evaluateObjective(obj.DesignVariables.x0)));
                    axis off;
                    
                case 'stress'
                    if isa(obj.Analysis, 'NonlinearAnalysis')
                        [~, stress, ~] = obj.Analysis.solveNonlinear(obj.Analysis.LoadVector);
                    else
                        stress = obj.Analysis.calculateStress();
                    end
                    surf(X, Y, stress.vm);
                    colormap('jet');
                    colorbar;
                    title('优化后von Mises应力');
            end
            
            % Save figure
            saveas(gcf, fullfile(obj.OutputFolder, 'OptimizationResults.png'));
            
            % Export optimal design
            obj.exportOptimalDesign();
        end
        
        function exportOptimalDesign(obj)
            % Export optimal design parameters and results
            filename = fullfile(obj.OutputFolder, 'optimal_design.txt');
            fid = fopen(filename, 'w');
            
            % Write header
            fprintf(fid, '优化设计结果\n');
            fprintf(fid, '=================\n\n');
            
            % Write objective type and value
            fprintf(fid, '目标函数类型: %s\n', obj.ObjectiveType.type);
            fprintf(fid, '最优目标函数值: %e\n\n', obj.evaluateObjective(obj.DesignVariables.x0));
            
            % Write design variables
            fprintf(fid, '设计变量:\n');
            idx = 1;
            for i = 1:length(obj.DesignVariables)
                var = obj.DesignVariables{i};
                fprintf(fid, '%s:\n', var.type);
                switch var.type
                    case {'thickness', 'radius'}
                        fprintf(fid, '  值: %e\n', obj.DesignVariables.x0(idx));
                        idx = idx + 1;
                    case 'shape'
                        n_params = length(var.control_points);
                        fprintf(fid, '  控制点参数:\n');
                        for j = 1:n_params
                            fprintf(fid, '    %d: %e\n', j, obj.DesignVariables.x0(idx+j-1));
                        end
                        idx = idx + n_params;
                end
            end
            
            % Write constraints
            if ~isempty(obj.Constraints)
                fprintf(fid, '\n约束条件:\n');
                [c, ceq] = obj.evaluateConstraints(obj.DesignVariables.x0);
                for i = 1:length(obj.Constraints)
                    fprintf(fid, '%s: %e\n', obj.Constraints{i}.type, c(i));
                end
            end
            
            fclose(fid);
        end
    end
end
