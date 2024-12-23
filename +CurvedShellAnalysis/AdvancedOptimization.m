classdef AdvancedOptimization < handle
    % ADVANCEDOPTIMIZATION Class for advanced optimization algorithms
    
    properties
        Problem  % Optimization problem definition
        Algorithm  % Algorithm type ('genetic', 'pso', 'surrogate', 'multi')
        Options  % Algorithm-specific options
        Results  % Optimization results
        OutputFolder = 'OptimizationResults'  % Output folder
    end
    
    methods
        function obj = AdvancedOptimization(problem, algorithm, options)
            % Constructor
            obj.Problem = problem;
            obj.Algorithm = algorithm;
            obj.Options = options;
            
            % Create output folder
            if ~exist(obj.OutputFolder, 'dir')
                mkdir(obj.OutputFolder);
            end
        end
        
        function [x_opt, f_opt, history] = optimize(obj)
            % Perform optimization using selected algorithm
            switch obj.Algorithm
                case 'genetic'
                    [x_opt, f_opt, history] = obj.geneticAlgorithm();
                case 'pso'
                    [x_opt, f_opt, history] = obj.particleSwarm();
                case 'surrogate'
                    [x_opt, f_opt, history] = obj.surrogateOptimization();
                case 'multi'
                    [x_opt, f_opt, history] = obj.multiObjective();
            end
            
            % Store results
            obj.Results = struct('x_opt', x_opt, 'f_opt', f_opt, 'history', history);
            
            % Plot results
            obj.plotResults();
        end
        
        function [x_opt, f_opt, history] = geneticAlgorithm(obj)
            % Genetic Algorithm implementation
            fprintf('运行遗传算法优化...\n');
            
            % Get problem parameters
            n_var = length(obj.Problem.lb);
            pop_size = obj.Options.pop_size;
            n_gen = obj.Options.n_generations;
            p_cross = obj.Options.p_crossover;
            p_mut = obj.Options.p_mutation;
            
            % Initialize population
            pop = obj.initializePopulation(pop_size, n_var);
            fit = zeros(pop_size, 1);
            
            % Initialize history
            history.best_fit = zeros(n_gen, 1);
            history.avg_fit = zeros(n_gen, 1);
            history.diversity = zeros(n_gen, 1);
            
            % Main loop
            for gen = 1:n_gen
                % Evaluate fitness
                parfor i = 1:pop_size
                    fit(i) = obj.evaluateFitness(pop(i,:));
                end
                
                % Store statistics
                [best_fit, best_idx] = min(fit);
                history.best_fit(gen) = best_fit;
                history.avg_fit(gen) = mean(fit);
                history.diversity(gen) = std(fit);
                
                % Selection
                parents = obj.tournamentSelection(pop, fit);
                
                % Crossover
                offspring = obj.crossover(parents, p_cross);
                
                % Mutation
                offspring = obj.mutation(offspring, p_mut);
                
                % Elitism
                pop = [pop(best_idx,:); offspring(1:end-1,:)];
            end
            
            % Get best solution
            [f_opt, best_idx] = min(fit);
            x_opt = pop(best_idx,:);
        end
        
        function [x_opt, f_opt, history] = particleSwarm(obj)
            % Particle Swarm Optimization implementation
            fprintf('运行粒子群优化...\n');
            
            % Get problem parameters
            n_var = length(obj.Problem.lb);
            swarm_size = obj.Options.swarm_size;
            n_iter = obj.Options.n_iterations;
            w = obj.Options.inertia;
            c1 = obj.Options.cognitive;
            c2 = obj.Options.social;
            
            % Initialize particles
            particles = obj.initializeParticles(swarm_size, n_var);
            velocities = zeros(swarm_size, n_var);
            p_best = particles;
            p_best_fit = inf(swarm_size, 1);
            g_best = zeros(1, n_var);
            g_best_fit = inf;
            
            % Initialize history
            history.best_fit = zeros(n_iter, 1);
            history.avg_fit = zeros(n_iter, 1);
            history.diversity = zeros(n_iter, 1);
            
            % Main loop
            for iter = 1:n_iter
                % Evaluate fitness
                parfor i = 1:swarm_size
                    fit = obj.evaluateFitness(particles(i,:));
                    
                    % Update personal best
                    if fit < p_best_fit(i)
                        p_best_fit(i) = fit;
                        p_best(i,:) = particles(i,:);
                        
                        % Update global best
                        if fit < g_best_fit
                            g_best_fit = fit;
                            g_best = particles(i,:);
                        end
                    end
                end
                
                % Store statistics
                history.best_fit(iter) = g_best_fit;
                history.avg_fit(iter) = mean(p_best_fit);
                history.diversity(iter) = std(p_best_fit);
                
                % Update velocities and positions
                r1 = rand(swarm_size, n_var);
                r2 = rand(swarm_size, n_var);
                
                velocities = w*velocities + ...
                    c1*r1.*(p_best - particles) + ...
                    c2*r2.*(repmat(g_best, swarm_size, 1) - particles);
                
                particles = particles + velocities;
                
                % Enforce bounds
                particles = max(min(particles, obj.Problem.ub), obj.Problem.lb);
            end
            
            % Get best solution
            x_opt = g_best;
            f_opt = g_best_fit;
        end
        
        function [x_opt, f_opt, history] = surrogateOptimization(obj)
            % Surrogate-based optimization implementation
            fprintf('运行代理模型优化...\n');
            
            % Get problem parameters
            n_var = length(obj.Problem.lb);
            n_initial = obj.Options.n_initial;
            n_iter = obj.Options.n_iterations;
            
            % Initialize sample points
            X = obj.latinHypercubeSampling(n_initial, n_var);
            y = zeros(n_initial, 1);
            
            % Evaluate initial points
            parfor i = 1:n_initial
                y(i) = obj.evaluateFitness(X(i,:));
            end
            
            % Initialize history
            history.best_fit = zeros(n_iter, 1);
            history.pred_error = zeros(n_iter, 1);
            
            % Main loop
            for iter = 1:n_iter
                % Train surrogate model
                model = obj.trainSurrogateModel(X, y);
                
                % Find next point to evaluate
                x_new = obj.findNextPoint(model, X);
                
                % Evaluate new point
                y_new = obj.evaluateFitness(x_new);
                
                % Update database
                X = [X; x_new];
                y = [y; y_new];
                
                % Store statistics
                [history.best_fit(iter), best_idx] = min(y);
                history.pred_error(iter) = obj.validateModel(model, X, y);
            end
            
            % Get best solution
            [f_opt, best_idx] = min(y);
            x_opt = X(best_idx,:);
        end
        
        function [x_opt, f_opt, history] = multiObjective(obj)
            % Multi-objective optimization using NSGA-II
            fprintf('运行多目标优化...\n');
            
            % Get problem parameters
            n_var = length(obj.Problem.lb);
            n_obj = length(obj.Problem.objectives);
            pop_size = obj.Options.pop_size;
            n_gen = obj.Options.n_generations;
            
            % Initialize population
            pop = obj.initializePopulation(pop_size, n_var);
            objectives = zeros(pop_size, n_obj);
            
            % Initialize history
            history.pareto_front = cell(n_gen, 1);
            history.diversity = zeros(n_gen, 1);
            
            % Main loop
            for gen = 1:n_gen
                % Evaluate objectives
                parfor i = 1:pop_size
                    objectives(i,:) = obj.evaluateObjectives(pop(i,:));
                end
                
                % Non-dominated sorting
                [ranks, crowding] = obj.nonDominatedSort(objectives);
                
                % Selection
                parents = obj.crowdedTournamentSelection(pop, ranks, crowding);
                
                % Crossover and mutation
                offspring = obj.geneticOperators(parents);
                
                % Evaluate offspring
                offspring_obj = zeros(size(offspring,1), n_obj);
                parfor i = 1:size(offspring,1)
                    offspring_obj(i,:) = obj.evaluateObjectives(offspring(i,:));
                end
                
                % Combine populations
                pop = [pop; offspring];
                objectives = [objectives; offspring_obj];
                
                % Select next generation
                [pop, objectives] = obj.environmentalSelection(pop, objectives);
                
                % Store statistics
                history.pareto_front{gen} = objectives(ranks==1,:);
                history.diversity(gen) = obj.calculateDiversity(objectives);
            end
            
            % Get Pareto optimal solutions
            [ranks, ~] = obj.nonDominatedSort(objectives);
            pareto_idx = ranks == 1;
            x_opt = pop(pareto_idx,:);
            f_opt = objectives(pareto_idx,:);
        end
        
        function pop = initializePopulation(obj, pop_size, n_var)
            % Initialize population with Latin Hypercube Sampling
            pop = obj.latinHypercubeSampling(pop_size, n_var);
        end
        
        function X = latinHypercubeSampling(obj, n_samples, n_var)
            % Generate Latin Hypercube samples
            X = zeros(n_samples, n_var);
            for i = 1:n_var
                X(:,i) = (randperm(n_samples) - 0.5)/n_samples;
            end
            
            % Scale to bounds
            X = X .* (obj.Problem.ub - obj.Problem.lb) + obj.Problem.lb;
        end
        
        function model = trainSurrogateModel(obj, X, y)
            % Train Gaussian Process surrogate model
            model.X = X;
            model.y = y;
            
            % Optimize hyperparameters
            theta0 = ones(1, size(X,2));
            model.theta = obj.optimizeHyperparameters(theta0, X, y);
        end
        
        function x_new = findNextPoint(obj, model, X)
            % Find next point using Expected Improvement
            n_var = size(X,2);
            
            % Define acquisition function
            acq_fun = @(x) obj.expectedImprovement(x, model);
            
            % Optimize acquisition function
            options = optimoptions('ga', 'Display', 'off');
            x_new = ga(@(x) -acq_fun(x), n_var, [], [], [], [], ...
                obj.Problem.lb, obj.Problem.ub, [], options);
        end
        
        function ei = expectedImprovement(obj, x, model)
            % Calculate Expected Improvement
            [mu, sigma] = obj.predictGP(x, model);
            y_min = min(model.y);
            
            % Handle numerical issues
            if sigma < 1e-6
                ei = 0;
                return
            end
            
            z = (y_min - mu)/sigma;
            ei = sigma * (z * normcdf(z) + normpdf(z));
        end
        
        function [ranks, crowding] = nonDominatedSort(obj, objectives)
            % Non-dominated sorting for multi-objective optimization
            n_pop = size(objectives, 1);
            n_obj = size(objectives, 2);
            
            % Initialize
            ranks = zeros(n_pop, 1);
            crowding = zeros(n_pop, 1);
            
            % Calculate domination
            domination = false(n_pop);
            for i = 1:n_pop
                for j = i+1:n_pop
                    if all(objectives(i,:) <= objectives(j,:)) && ...
                            any(objectives(i,:) < objectives(j,:))
                        domination(i,j) = true;
                    elseif all(objectives(j,:) <= objectives(i,:)) && ...
                            any(objectives(j,:) < objectives(i,:))
                        domination(j,i) = true;
                    end
                end
            end
            
            % Assign ranks
            rank = 1;
            while any(ranks == 0)
                non_dominated = find(ranks == 0 & ~any(domination(:,ranks == 0), 1))';
                ranks(non_dominated) = rank;
                rank = rank + 1;
            end
            
            % Calculate crowding distance
            for r = 1:max(ranks)
                idx = ranks == r;
                n_idx = sum(idx);
                
                if n_idx <= 2
                    crowding(idx) = inf;
                    continue;
                end
                
                % Sort objectives
                sorted_obj = sort(objectives(idx,:));
                
                % Calculate distances
                for j = 1:n_obj
                    [~, order] = sort(objectives(idx,j));
                    crowding(idx) = crowding(idx) + ...
                        (sorted_obj(end,j) - sorted_obj(1,j)) .* ...
                        [inf; diff(objectives(idx(order),j)); inf];
                end
            end
        end
        
        function diversity = calculateDiversity(obj, objectives)
            % Calculate diversity metric for Pareto front
            % Using crowding distance
            n_obj = size(objectives, 2);
            
            % Normalize objectives
            obj_min = min(objectives);
            obj_max = max(objectives);
            obj_norm = (objectives - obj_min) ./ (obj_max - obj_min);
            
            % Calculate average distance to nearest neighbor
            distances = zeros(size(objectives,1), 1);
            for i = 1:size(objectives,1)
                dist = sqrt(sum((obj_norm - obj_norm(i,:)).^2, 2));
                dist(i) = inf;
                distances(i) = min(dist);
            end
            
            diversity = mean(distances);
        end
        
        function plotResults(obj)
            % Plot optimization results
            switch obj.Algorithm
                case {'genetic', 'pso', 'surrogate'}
                    obj.plotSingleObjective();
                case 'multi'
                    obj.plotMultiObjective();
            end
        end
        
        function plotSingleObjective(obj)
            % Plot single-objective optimization results
            figure('Position', [100 100 800 600]);
            
            % Convergence history
            subplot(2,1,1);
            plot(obj.Results.history.best_fit, 'b-', 'LineWidth', 2);
            hold on;
            plot(obj.Results.history.avg_fit, 'r--', 'LineWidth', 1);
            grid on;
            title('优化收敛历史');
            xlabel('迭代次数');
            ylabel('目标函数值');
            legend('最优值', '平均值');
            
            % Population diversity
            subplot(2,1,2);
            plot(obj.Results.history.diversity, 'g-', 'LineWidth', 2);
            grid on;
            title('种群多样性');
            xlabel('迭代次数');
            ylabel('多样性指标');
            
            % Save figure
            saveas(gcf, fullfile(obj.OutputFolder, 'SingleObjectiveResults.png'));
        end
        
        function plotMultiObjective(obj)
            % Plot multi-objective optimization results
            figure('Position', [100 100 800 600]);
            
            % Final Pareto front
            subplot(2,1,1);
            scatter(obj.Results.f_opt(:,1), obj.Results.f_opt(:,2), 'filled');
            grid on;
            title('Pareto前沿');
            xlabel('目标函数1');
            ylabel('目标函数2');
            
            % Diversity history
            subplot(2,1,2);
            plot(obj.Results.history.diversity, 'g-', 'LineWidth', 2);
            grid on;
            title('Pareto前沿多样性');
            xlabel('迭代次数');
            ylabel('多样性指标');
            
            % Save figure
            saveas(gcf, fullfile(obj.OutputFolder, 'MultiObjectiveResults.png'));
        end
    end
end
