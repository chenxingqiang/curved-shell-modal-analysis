classdef ParallelManager < handle
    % PARALLELMANAGER Class for parallel computing management
    
    properties
        NumWorkers  % Number of parallel workers
        BatchSize  % Batch size for parallel processing
        MemoryLimit  % Memory limit per worker
        GPUDevices  % Available GPU devices
        ClusterProfile  % Parallel cluster profile
        TaskQueue  % Queue of parallel tasks
        Results  % Results from parallel execution
    end
    
    methods
        function obj = ParallelManager(num_workers, batch_size, memory_limit)
            % Constructor
            obj.NumWorkers = num_workers;
            obj.BatchSize = batch_size;
            obj.MemoryLimit = memory_limit;
            obj.TaskQueue = {};
            obj.Results = {};
            
            % Initialize parallel pool
            obj.initializeParallel();
        end
        
        function initializeParallel(obj)
            % Initialize parallel computing environment
            fprintf('初始化并行计算环境...\n');
            
            % Check if parallel pool exists
            if isempty(gcp('nocreate'))
                % Create parallel pool
                parpool('local', obj.NumWorkers);
            end
            
            % Check GPU availability
            if gpuDeviceCount > 0
                obj.GPUDevices = gpuDevice;
                fprintf('找到 %d 个GPU设备\n', length(obj.GPUDevices));
            else
                obj.GPUDevices = [];
                fprintf('未找到GPU设备\n');
            end
            
            % Set memory limit
            if ~isempty(obj.MemoryLimit)
                parallel.pool.Constant('MemLimit', obj.MemoryLimit);
            end
        end
        
        function addTask(obj, task_type, params)
            % Add task to queue
            task = struct('type', task_type, 'params', params);
            obj.TaskQueue{end+1} = task;
        end
        
        function results = executeTasks(obj)
            % Execute tasks in parallel
            fprintf('执行并行任务...\n');
            
            % Group tasks into batches
            n_tasks = length(obj.TaskQueue);
            n_batches = ceil(n_tasks/obj.BatchSize);
            
            % Process batches
            for i = 1:n_batches
                start_idx = (i-1)*obj.BatchSize + 1;
                end_idx = min(i*obj.BatchSize, n_tasks);
                batch_tasks = obj.TaskQueue(start_idx:end_idx);
                
                % Execute batch in parallel
                batch_results = obj.executeBatch(batch_tasks);
                
                % Store results
                obj.Results = [obj.Results batch_results];
            end
            
            % Return results
            results = obj.Results;
            
            % Clear task queue
            obj.TaskQueue = {};
        end
        
        function results = executeBatch(obj, batch_tasks)
            % Execute a batch of tasks in parallel
            n_tasks = length(batch_tasks);
            results = cell(1, n_tasks);
            
            % Determine execution method
            if ~isempty(obj.GPUDevices) && ...
                    all(cellfun(@(x) strcmp(x.type, 'gpu_compatible'), batch_tasks))
                % GPU parallel execution
                results = obj.executeGPUBatch(batch_tasks);
            else
                % CPU parallel execution
                parfor i = 1:n_tasks
                    results{i} = obj.executeTask(batch_tasks{i});
                end
            end
        end
        
        function result = executeTask(obj, task)
            % Execute a single task
            switch task.type
                case 'modal_analysis'
                    result = obj.executeModalAnalysis(task.params);
                case 'nonlinear_analysis'
                    result = obj.executeNonlinearAnalysis(task.params);
                case 'optimization'
                    result = obj.executeOptimization(task.params);
                case 'material_model'
                    result = obj.executeMaterialModel(task.params);
            end
        end
        
        function results = executeGPUBatch(obj, batch_tasks)
            % Execute tasks on GPU
            n_tasks = length(batch_tasks);
            n_gpus = length(obj.GPUDevices);
            results = cell(1, n_tasks);
            
            % Distribute tasks across GPUs
            spmd
                % Get GPU device for this worker
                gpu_idx = mod(labindex-1, n_gpus) + 1;
                gpuDevice(gpu_idx);
                
                % Get tasks for this worker
                worker_tasks = batch_tasks(labindex:numlabs:end);
                
                % Execute tasks
                worker_results = cell(size(worker_tasks));
                for i = 1:length(worker_tasks)
                    worker_results{i} = obj.executeGPUTask(worker_tasks{i});
                end
            end
            
            % Combine results
            results = [worker_results{:}];
        end
        
        function result = executeGPUTask(obj, task)
            % Execute a single task on GPU
            switch task.type
                case 'modal_analysis'
                    result = obj.executeModalAnalysisGPU(task.params);
                case 'nonlinear_analysis'
                    result = obj.executeNonlinearAnalysisGPU(task.params);
                case 'optimization'
                    result = obj.executeOptimizationGPU(task.params);
            end
        end
        
        function result = executeModalAnalysis(obj, params)
            % Execute modal analysis task
            surface = params.surface;
            n_modes = params.n_modes;
            
            % Create analysis object
            analysis = CurvedShellAnalysis.ModalAnalysis(surface, n_modes);
            
            % Run analysis
            analysis.analyze();
            
            % Return results
            result = struct('frequencies', analysis.NaturalFrequencies, ...
                          'mode_shapes', analysis.ModeShapes);
        end
        
        function result = executeModalAnalysisGPU(obj, params)
            % Execute modal analysis task on GPU
            surface = params.surface;
            n_modes = params.n_modes;
            
            % Transfer matrices to GPU
            K = gpuArray(surface.StiffnessMatrix);
            M = gpuArray(surface.MassMatrix);
            
            % Solve eigenvalue problem on GPU
            [V, D] = eigs(K, M, n_modes, 'smallestabs');
            
            % Transfer results back to CPU
            frequencies = sqrt(diag(gather(D)))/(2*pi);
            mode_shapes = gather(V);
            
            % Return results
            result = struct('frequencies', frequencies, ...
                          'mode_shapes', mode_shapes);
        end
        
        function result = executeNonlinearAnalysis(obj, params)
            % Execute nonlinear analysis task
            surface = params.surface;
            load_steps = params.load_steps;
            
            % Create analysis object
            analysis = CurvedShellAnalysis.NonlinearAnalysis(surface);
            
            % Run analysis
            [U, stress, strain] = analysis.solveNonlinear(load_steps);
            
            % Return results
            result = struct('displacement', U, ...
                          'stress', stress, ...
                          'strain', strain);
        end
        
        function result = executeNonlinearAnalysisGPU(obj, params)
            % Execute nonlinear analysis task on GPU
            surface = params.surface;
            load_steps = params.load_steps;
            
            % Transfer matrices to GPU
            K = gpuArray(surface.StiffnessMatrix);
            F = gpuArray(load_steps);
            
            % Initialize
            n = size(K,1);
            U = gpuArray(zeros(n, 1));
            
            % Newton-Raphson iteration on GPU
            for step = 1:length(load_steps)
                R = F(:,step) - K*U;
                dU = K\R;
                U = U + dU;
            end
            
            % Calculate stress and strain on GPU
            [stress, strain] = surface.calculateStressStrainGPU(gather(U));
            
            % Return results
            result = struct('displacement', gather(U), ...
                          'stress', stress, ...
                          'strain', strain);
        end
        
        function result = executeOptimization(obj, params)
            % Execute optimization task
            problem = params.problem;
            algorithm = params.algorithm;
            options = params.options;
            
            % Create optimization object
            opt = CurvedShellAnalysis.AdvancedOptimization(problem, algorithm, options);
            
            % Run optimization
            [x_opt, f_opt, history] = opt.optimize();
            
            % Return results
            result = struct('x_opt', x_opt, ...
                          'f_opt', f_opt, ...
                          'history', history);
        end
        
        function result = executeOptimizationGPU(obj, params)
            % Execute optimization task on GPU
            problem = params.problem;
            algorithm = params.algorithm;
            options = params.options;
            
            % Transfer data to GPU
            problem.objective = @(x) gather(problem.objective(gpuArray(x)));
            if isfield(problem, 'constraint')
                problem.constraint = @(x) gather(problem.constraint(gpuArray(x)));
            end
            
            % Create optimization object
            opt = CurvedShellAnalysis.AdvancedOptimization(problem, algorithm, options);
            
            % Run optimization
            [x_opt, f_opt, history] = opt.optimize();
            
            % Return results
            result = struct('x_opt', x_opt, ...
                          'f_opt', f_opt, ...
                          'history', history);
        end
        
        function result = executeMaterialModel(obj, params)
            % Execute material model task
            material_type = params.type;
            properties = params.properties;
            strain = params.strain;
            
            % Create material object
            material = CurvedShellAnalysis.AdvancedMaterial(material_type, properties);
            
            % Calculate response
            [stress, D, state] = material.calculateStress(strain);
            
            % Return results
            result = struct('stress', stress, ...
                          'tangent', D, ...
                          'state', state);
        end
    end
end
