classdef BucklingAnalysis < handle
    % BUCKLINGANALYSIS Class for linear and nonlinear buckling analysis
    
    properties
        Surface  % Surface object
        ModalAnalysis  % Reference to modal analysis object
        LoadType  % Type of loading ('compression', 'shear', 'pressure')
        LoadMagnitude  % Load magnitude
        BucklingModes  % Buckling mode shapes
        BucklingLoads  % Critical buckling loads
        OutputFolder = 'BucklingResults'  % Output folder
    end
    
    methods
        function obj = BucklingAnalysis(modalAnalysis)
            % Constructor
            obj.ModalAnalysis = modalAnalysis;
            obj.Surface = modalAnalysis.Surface;
            
            % Create output folder
            if ~exist(obj.OutputFolder, 'dir')
                mkdir(obj.OutputFolder);
            end
        end
        
        function setLoading(obj, type, magnitude)
            % Set loading conditions
            obj.LoadType = type;
            obj.LoadMagnitude = magnitude;
        end
        
        function [lambda, modes] = linearBuckling(obj, num_modes)
            % Perform linear buckling analysis
            fprintf('执行线性屈曲分析...\n');
            
            % Get stiffness matrix
            K = obj.ModalAnalysis.StiffnessMatrix;
            
            % Calculate geometric stiffness matrix
            Kg = obj.calculateGeometricStiffness();
            
            % Solve generalized eigenvalue problem
            [V, D] = eigs(K, Kg, num_modes, 'smallestabs');
            
            % Extract buckling loads and modes
            obj.BucklingLoads = diag(D);
            obj.BucklingModes = V;
            
            % Sort by magnitude
            [obj.BucklingLoads, idx] = sort(obj.BucklingLoads);
            obj.BucklingModes = obj.BucklingModes(:,idx);
            
            % Return results
            lambda = obj.BucklingLoads;
            modes = obj.BucklingModes;
            
            % Visualize results
            obj.plotBucklingModes();
        end
        
        function Kg = calculateGeometricStiffness(obj)
            % Calculate geometric stiffness matrix
            [X, Y, Z] = obj.Surface.generateMesh();
            nx = obj.Surface.nx;
            ny = obj.Surface.ny;
            ndof = obj.ModalAnalysis.DOFsPerNode;
            N = nx * ny * ndof;
            
            % Initialize geometric stiffness matrix
            Kg = sparse(N, N);
            
            % Assembly loop
            for i = 1:nx-1
                for j = 1:ny-1
                    % Get element nodes
                    nodes = [j+(i-1)*ny, j+i*ny, j+1+i*ny, j+1+(i-1)*ny];
                    
                    % Calculate element geometric stiffness
                    Kg_e = obj.calculateElementGeometricStiffness(X(nodes), Y(nodes), Z(nodes));
                    
                    % Assemble into global matrix
                    dofs = [];
                    for n = nodes
                        dofs = [dofs, (n-1)*ndof+1:n*ndof];
                    end
                    Kg(dofs,dofs) = Kg(dofs,dofs) + Kg_e;
                end
            end
        end
        
        function Kg_e = calculateElementGeometricStiffness(obj, xe, ye, ze)
            % Calculate element geometric stiffness matrix
            % This is a simplified implementation
            
            % Element size
            dx = abs(xe(2) - xe(1));
            dy = abs(ye(4) - ye(1));
            
            % Pre-stress state based on load type
            switch lower(obj.LoadType)
                case 'compression'
                    sigma_x = -obj.LoadMagnitude;
                    sigma_y = 0;
                    tau_xy = 0;
                case 'shear'
                    sigma_x = 0;
                    sigma_y = 0;
                    tau_xy = obj.LoadMagnitude;
                case 'pressure'
                    sigma_x = -obj.LoadMagnitude/2;
                    sigma_y = -obj.LoadMagnitude/2;
                    tau_xy = 0;
            end
            
            % Basic geometric stiffness matrix for a shell element
            Kg_e = zeros(20, 20);  % 5 DOF per node, 4 nodes
            
            % Simplified geometric stiffness terms
            k = [sigma_x*dy/dx   tau_xy
                 tau_xy         sigma_y*dx/dy];
            
            % Fill the geometric stiffness matrix
            for i = 1:4
                for j = 1:4
                    idx_i = (i-1)*5 + (1:2);
                    idx_j = (j-1)*5 + (1:2);
                    Kg_e(idx_i,idx_j) = k;
                end
            end
        end
        
        function plotBucklingModes(obj)
            % Plot buckling mode shapes
            [X, Y, Z] = obj.Surface.generateMesh();
            nx = obj.Surface.nx;
            ny = obj.Surface.ny;
            
            for mode = 1:size(obj.BucklingModes,2)
                figure('Position', [100 100 800 600]);
                
                % Extract modal displacements
                w = zeros(ny, nx);
                for i = 1:nx
                    for j = 1:ny
                        n = j + (i-1)*ny;
                        w(j,i) = obj.BucklingModes((n-1)*5+3, mode);
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
                title(sprintf('屈曲模态 %d (λ = %.2f)', mode, obj.BucklingLoads(mode)));
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
                plot(X(ceil(ny/2),:), Zdef(ceil(ny/2),:), 'b-', 'LineWidth', 2);
                hold on;
                plot(X(ceil(ny/2),:), Z(ceil(ny/2),:), 'k--', 'LineWidth', 1);
                title('中心线变形');
                xlabel('X (m)');
                ylabel('Z (m)');
                legend('变形', '原始形状');
                axis equal;
                grid on;
                
                % Save figure
                saveas(gcf, fullfile(obj.OutputFolder, sprintf('BucklingMode%d.png', mode)));
            end
            
            % Export to STL
            obj.exportToSTL(X, Y, Z, w);
        end
        
        function exportToSTL(obj, X, Y, Z, w)
            % Export deformed shape to STL format
            Zdef = Z + w*0.1*max(abs(Z(:)));
            
            % Create triangulation
            [nx, ny] = size(X);
            tri = [];
            for i = 1:nx-1
                for j = 1:ny-1
                    % First triangle
                    tri = [tri; [i+(j-1)*nx, (i+1)+(j-1)*nx, i+j*nx]];
                    % Second triangle
                    tri = [tri; [(i+1)+(j-1)*nx, (i+1)+j*nx, i+j*nx]];
                end
            end
            
            % Combine coordinates
            vertices = [X(:) Y(:) Zdef(:)];
            
            % Write STL file
            filename = fullfile(obj.OutputFolder, 'buckling_mode.stl');
            obj.writeSTL(filename, vertices, tri);
        end
        
        function writeSTL(~, filename, vertices, triangles)
            % Write STL file
            fid = fopen(filename, 'w');
            
            % Write header
            fprintf(fid, 'solid BucklingMode\n');
            
            % Write triangles
            for i = 1:size(triangles, 1)
                % Get vertex coordinates
                v1 = vertices(triangles(i,1),:);
                v2 = vertices(triangles(i,2),:);
                v3 = vertices(triangles(i,3),:);
                
                % Calculate normal vector
                normal = cross(v2-v1, v3-v1);
                normal = normal/norm(normal);
                
                % Write facet
                fprintf(fid, 'facet normal %e %e %e\n', normal);
                fprintf(fid, '  outer loop\n');
                fprintf(fid, '    vertex %e %e %e\n', v1);
                fprintf(fid, '    vertex %e %e %e\n', v2);
                fprintf(fid, '    vertex %e %e %e\n', v3);
                fprintf(fid, '  endloop\n');
                fprintf(fid, 'endfacet\n');
            end
            
            % Write footer
            fprintf(fid, 'endsolid BucklingMode\n');
            fclose(fid);
        end
    end
end
