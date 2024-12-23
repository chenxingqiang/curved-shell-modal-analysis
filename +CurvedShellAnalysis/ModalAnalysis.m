classdef ModalAnalysis < handle
    % MODALANALYSIS Class for performing modal analysis on curved surfaces
    
    properties
        Surface  % Surface object
        NumModes = 6  % Number of modes to compute
        OutputFolder = 'ModalResults'  % Output folder for results
        StiffnessMatrix  % Global stiffness matrix
        MassMatrix  % Global mass matrix
        DampingMatrix  % Global damping matrix
        Frequencies  % Natural frequencies
        ModesMatrix  % Mode shapes
        DOFsPerNode = 5  % Degrees of freedom per node
    end
    
    methods
        function obj = ModalAnalysis(surface, numModes)
            % Constructor
            obj.Surface = surface;
            if nargin > 1
                obj.NumModes = numModes;
            end
            
            % Create output folder if it doesn't exist
            if ~exist(obj.OutputFolder, 'dir')
                mkdir(obj.OutputFolder);
            end
        end
        
        function analyze(obj)
            % Perform complete modal analysis
            fprintf('开始模态分析...\n');
            
            % Generate mesh
            [X, Y, Z] = obj.Surface.generateMesh();
            
            % Assemble matrices
            obj.assembleMatrices(X, Y, Z);
            
            % Apply boundary conditions and solve
            obj.solveEigenvalueProblem();
            
            % Save results
            obj.saveResults();
            
            % Visualize results
            obj.visualizeResults(X, Y, Z);
            
            fprintf('分析完成！结果已保存到 %s\n', obj.OutputFolder);
        end
        
        function assembleMatrices(obj, X, Y, Z)
            % Assemble global matrices
            fprintf('正在组装全局矩阵...\n');
            
            nx = obj.Surface.nx;
            ny = obj.Surface.ny;
            N = nx * ny * obj.DOFsPerNode;
            
            % Initialize sparse matrices
            obj.StiffnessMatrix = sparse(N, N);
            obj.MassMatrix = sparse(N, N);
            
            h = waitbar(0, '组装单元矩阵...');
            
            % Assembly loop
            for i = 1:nx-1
                for j = 1:ny-1
                    % Update progress bar
                    waitbar(((i-1)*(ny-1) + j)/((nx-1)*(ny-1)), h);
                    
                    % Get element nodes
                    nodes = [j+(i-1)*ny, j+i*ny, j+1+i*ny, j+1+(i-1)*ny];
                    
                    % Get element coordinates
                    xe = X(nodes);
                    ye = Y(nodes);
                    ze = Z(nodes);
                    
                    % Calculate element matrices
                    [ke, me] = obj.calculateElementMatrices(xe, ye, ze);
                    
                    % Assemble into global matrices
                    dofs = [];
                    for n = nodes
                        dofs = [dofs, (n-1)*obj.DOFsPerNode+1:n*obj.DOFsPerNode];
                    end
                    
                    obj.StiffnessMatrix(dofs,dofs) = obj.StiffnessMatrix(dofs,dofs) + ke;
                    obj.MassMatrix(dofs,dofs) = obj.MassMatrix(dofs,dofs) + me;
                end
            end
            
            close(h);
        end
        
        function [ke, me] = calculateElementMatrices(obj, xe, ye, ze)
            % Calculate element stiffness and mass matrices
            % This is a simplified implementation
            
            % Element size
            dx = abs(xe(2) - xe(1));
            dy = abs(ye(4) - ye(1));
            
            % Material properties
            E = obj.Surface.E;
            nu = obj.Surface.nu;
            t = obj.Surface.t;
            rho = obj.Surface.rho;
            
            % Membrane stiffness
            ke_membrane = E*t/(1-nu^2) * dx*dy/4 * [
                4  2  0  2  1  0  1  0  0  2  1  0
                2  4  0  1  2  0  0  1  0  1  2  0
                0  0  1  0  0  0  0  0  0  0  0  0
                2  1  0  4  2  0  2  1  0  1  0  0
                1  2  0  2  4  0  1  2  0  0  1  0
                0  0  0  0  0  1  0  0  0  0  0  0
                1  0  0  2  1  0  4  2  0  2  1  0
                0  1  0  1  2  0  2  4  0  1  2  0
                0  0  0  0  0  0  0  0  1  0  0  0
                2  1  0  1  0  0  2  1  0  4  2  0
                1  2  0  0  1  0  1  2  0  2  4  0
                0  0  0  0  0  0  0  0  0  0  0  1
            ];
            
            % Bending stiffness
            D = E*t^3/(12*(1-nu^2));
            ke_bending = D*dx*dy/4 * [
                4  1  1  0  1  0  0  0
                1  4  0  1  0  1  0  0
                1  0  4  1  0  0  1  0
                0  1  1  4  0  0  0  1
                1  0  0  0  4  1  1  0
                0  1  0  0  1  4  0  1
                0  0  1  0  1  0  4  1
                0  0  0  1  0  1  1  4
            ];
            
            % Combine membrane and bending
            ke = zeros(20, 20);
            ke(1:12, 1:12) = ke_membrane;
            ke(13:20, 13:20) = ke_bending;
            
            % Mass matrix
            me = zeros(20, 20);
            
            % Translational mass
            me(1:12, 1:12) = rho*t*dx*dy/36 * [
                4  0  0  2  0  0  1  0  0  2  0  0
                0  4  0  0  2  0  0  1  0  0  2  0
                0  0  4  0  0  2  0  0  1  0  0  2
                2  0  0  4  0  0  2  0  0  1  0  0
                0  2  0  0  4  0  0  2  0  0  1  0
                0  0  2  0  0  4  0  0  2  0  0  1
                1  0  0  2  0  0  4  0  0  2  0  0
                0  1  0  0  2  0  0  4  0  0  2  0
                0  0  1  0  0  2  0  0  4  0  0  2
                2  0  0  1  0  0  2  0  0  4  0  0
                0  2  0  0  1  0  0  2  0  0  4  0
                0  0  2  0  0  1  0  0  2  0  0  4
            ];
            
            % Rotational mass
            me(13:20, 13:20) = rho*t^3*dx*dy/420 * [
                16  0 -8   0  4  0 -8   0
                 0 16  0  -8  0  4  0  -8
                -8  0 16   0 -8  0  4   0
                 0 -8  0  16  0 -8  0   4
                 4  0 -8   0 16  0 -8   0
                 0  4  0  -8  0 16  0  -8
                -8  0  4   0 -8  0 16   0
                 0 -8  0   4  0 -8  0  16
            ];
        end
        
        function solveEigenvalueProblem(obj)
            % Solve the eigenvalue problem
            fprintf('应用边界条件...\n');
            
            % Apply boundary conditions (fixed ends)
            nx = obj.Surface.nx;
            ny = obj.Surface.ny;
            fixed_dofs = [];
            for i = [1, nx]
                for j = 1:ny
                    n = j + (i-1)*ny;
                    fixed_dofs = [fixed_dofs, (n-1)*obj.DOFsPerNode+1:n*obj.DOFsPerNode];
                end
            end
            free_dofs = setdiff(1:size(obj.StiffnessMatrix,1), fixed_dofs);
            
            fprintf('求解特征值问题...\n');
            [V, D] = eig(full(obj.StiffnessMatrix(free_dofs,free_dofs)), ...
                        full(obj.MassMatrix(free_dofs,free_dofs)));
            
            % Extract frequencies and mode shapes
            [obj.Frequencies, idx] = sort(sqrt(real(diag(D)))/(2*pi));
            obj.Frequencies = obj.Frequencies(1:obj.NumModes);
            
            % Build complete mode shapes
            V = V(:,idx);
            obj.ModesMatrix = zeros(size(obj.StiffnessMatrix,1), obj.NumModes);
            for i = 1:obj.NumModes
                mode_free = V(:,i);
                mode_full = zeros(size(obj.StiffnessMatrix,1), 1);
                mode_full(free_dofs) = mode_free;
                obj.ModesMatrix(:,i) = mode_full;
            end
            
            % Calculate damping matrix
            zeta = obj.Surface.zeta;
            alpha = 2*zeta*sqrt(obj.Frequencies(1)*obj.Frequencies(end));
            beta = 2*zeta/(sqrt(obj.Frequencies(1)*obj.Frequencies(end)));
            obj.DampingMatrix = alpha*obj.MassMatrix + beta*obj.StiffnessMatrix;
        end
        
        function saveResults(obj)
            % Save analysis results to file
            fid = fopen(fullfile(obj.OutputFolder, 'NaturalFrequencies.txt'), 'w');
            
            % Write header
            fprintf(fid, '曲面板模态分析结果\n\n');
            
            % Write geometry parameters
            fprintf(fid, '几何参数:\n');
            fprintf(fid, '长度: %.3f m\n', obj.Surface.L);
            fprintf(fid, '宽度: %.3f m\n', obj.Surface.W);
            fprintf(fid, '厚度: %.3f m\n', obj.Surface.t);
            
            % Write surface-specific parameters
            if isa(obj.Surface, 'CurvedShellAnalysis.SphericalSurface')
                fprintf(fid, '曲率半径: %.3f m\n', obj.Surface.R);
            elseif isa(obj.Surface, 'CurvedShellAnalysis.EllipsoidalSurface')
                fprintf(fid, 'X方向曲率半径: %.3f m\n', obj.Surface.Rx);
                fprintf(fid, 'Y方向曲率半径: %.3f m\n', obj.Surface.Ry);
                fprintf(fid, 'Z方向曲率半径: %.3f m\n', obj.Surface.Rz);
            end
            
            % Write material properties
            fprintf(fid, '\n材料属性:\n');
            fprintf(fid, '弹性模量: %.2e Pa\n', obj.Surface.E);
            fprintf(fid, '密度: %.1f kg/m³\n', obj.Surface.rho);
            fprintf(fid, '泊松比: %.2f\n', obj.Surface.nu);
            fprintf(fid, '阻尼比: %.3f\n\n', obj.Surface.zeta);
            
            % Write frequencies
            fprintf(fid, '前%d阶固有频率:\n', obj.NumModes);
            for i = 1:obj.NumModes
                fprintf(fid, '第%d阶: %.2f Hz\n', i, obj.Frequencies(i));
            end
            
            % Write damping parameters
            alpha = 2*obj.Surface.zeta*sqrt(obj.Frequencies(1)*obj.Frequencies(end));
            beta = 2*obj.Surface.zeta/(sqrt(obj.Frequencies(1)*obj.Frequencies(end)));
            fprintf(fid, '\n阻尼参数:\n');
            fprintf(fid, 'α (质量比例系数): %.6e\n', alpha);
            fprintf(fid, 'β (刚度比例系数): %.6e\n', beta);
            
            fclose(fid);
        end
        
        function visualizeResults(obj, X, Y, Z)
            % Visualize mode shapes and frequency response
            fprintf('绘制模态振型...\n');
            
            % Plot each mode shape
            for mode = 1:obj.NumModes
                figure('Position', [100 100 800 600]);
                
                % Extract modal displacements
                w = zeros(obj.Surface.ny, obj.Surface.nx);
                for i = 1:obj.Surface.nx
                    for j = 1:obj.Surface.ny
                        n = j + (i-1)*obj.Surface.ny;
                        w(j,i) = obj.ModesMatrix((n-1)*obj.DOFsPerNode+3, mode);
                    end
                end
                
                % Normalize displacement
                w = w/max(abs(w(:)));
                
                % Create deformed surface
                Zdef = Z + w*0.1*max(abs(Z(:)));
                
                % 3D surface plot
                subplot(2,2,[1,2]);
                surf(X, Y, Zdef);
                colormap('jet');
                shading interp;
                colorbar;
                title(sprintf('第%d阶模态 (f = %.2f Hz)', mode, obj.Frequencies(mode)));
                xlabel('X (m)');
                ylabel('Y (m)');
                zlabel('Z (m)');
                axis equal;
                view(45, 30);
                
                % Top view
                subplot(2,2,3);
                contourf(X, Y, w, 20);
                colormap('jet');
                colorbar;
                title('俯视图');
                xlabel('X (m)');
                ylabel('Y (m)');
                axis equal;
                
                % Side view
                subplot(2,2,4);
                plot(X(ceil(obj.Surface.ny/2),:), Zdef(ceil(obj.Surface.ny/2),:), ...
                     'b-', 'LineWidth', 2);
                hold on;
                plot(X(ceil(obj.Surface.ny/2),:), Z(ceil(obj.Surface.ny/2),:), ...
                     'k--', 'LineWidth', 1);
                title('中心线变形');
                xlabel('X (m)');
                ylabel('Z (m)');
                legend('变形', '原始形状');
                axis equal;
                grid on;
                
                % Save figure
                saveas(gcf, fullfile(obj.OutputFolder, sprintf('Mode%d.png', mode)));
            end
            
            % Calculate and plot frequency response function
            fprintf('计算频率响应函数...\n');
            f = logspace(0, 4, 1000);
            w = 2*pi*f;
            H = zeros(length(f), 1);
            
            % Select observation point (panel center)
            obs_node = ceil(obj.Surface.nx*obj.Surface.ny/2);
            obs_dof = (obs_node-1)*obj.DOFsPerNode + 3;
            
            % Calculate FRF
            for i = 1:length(f)
                Z = -w(i)^2*obj.MassMatrix + 1i*w(i)*obj.DampingMatrix + obj.StiffnessMatrix;
                F = zeros(size(Z,1), 1);
                F(obs_dof) = 1;
                U = Z\F;
                H(i) = abs(U(obs_dof));
            end
            
            % Plot FRF
            figure('Position', [100 100 800 400]);
            semilogx(f, 20*log10(H), 'LineWidth', 2);
            grid on;
            title('频率响应函数');
            xlabel('频率 (Hz)');
            ylabel('幅值 (dB)');
            saveas(gcf, fullfile(obj.OutputFolder, 'FRF.png'));
        end
    end
end
