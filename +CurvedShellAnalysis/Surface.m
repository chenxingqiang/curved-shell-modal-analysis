classdef Surface < handle
    % Surface Base class for all curved surfaces
    
    properties
        L           % Length
        W           % Width
        t           % Thickness
        E           % Young's modulus
        nu          % Poisson's ratio
        rho         % Density
        mesh        % Mesh object
        material    % Material properties
        nx = 20     % Number of elements in x direction
        ny = 20     % Number of elements in y direction
    end
    
    methods
        function obj = Surface(params)
            % Constructor for Surface class
            if nargin > 0
                if isfield(params, 'L')
                    obj.L = params.L;
                end
                if isfield(params, 'W')
                    obj.W = params.W;
                end
                if isfield(params, 't')
                    obj.t = params.t;
                end
                if isfield(params, 'E')
                    obj.E = params.E;
                end
                if isfield(params, 'nu')
                    obj.nu = params.nu;
                end
                if isfield(params, 'rho')
                    obj.rho = params.rho;
                end
                if isfield(params, 'nx')
                    obj.nx = params.nx;
                end
                if isfield(params, 'ny')
                    obj.ny = params.ny;
                end
            end
            obj.material = struct();
            obj.mesh = struct('X', [], 'Y', [], 'Z', []);
            obj.validateParams();
        end
        
        function setMaterial(obj, material)
            % Set material properties
            obj.material = material;
            if isfield(material, 'E')
                obj.E = material.E;
            end
            if isfield(material, 'nu')
                obj.nu = material.nu;
            end
            if isfield(material, 'rho')
                obj.rho = material.rho;
            end
        end
        
        function mesh = getMesh(obj)
            % Get mesh object
            mesh = obj.mesh;
        end
        
        function validateParams(obj)
            % Validate common parameters
            if ~isempty(obj.t)
                assert(obj.t > 0, 'Thickness must be positive');
            end
            if ~isempty(obj.E)
                assert(obj.E > 0, 'Young''s modulus must be positive');
            end
            if ~isempty(obj.nu)
                assert(obj.nu > -1 && obj.nu < 0.5, 'Poisson''s ratio must be between -1 and 0.5');
            end
            if ~isempty(obj.rho)
                assert(obj.rho > 0, 'Density must be positive');
            end
            if ~isempty(obj.nx)
                assert(obj.nx > 0, 'Number of elements in x direction must be positive');
            end
            if ~isempty(obj.ny)
                assert(obj.ny > 0, 'Number of elements in y direction must be positive');
            end
        end
        
        function plotScalarField(obj, field, title_str)
            % Plot a scalar field on the surface
            % field: 2D array matching mesh dimensions
            % title_str: Optional title string
            
            if nargin < 3
                title_str = 'Scalar Field';
            end
            
            % Check field dimensions
            [ny, nx] = size(obj.mesh.X);
            assert(all(size(field) == [ny, nx]), ...
                   'Field dimensions must match mesh dimensions');
            
            % Create surface plot with scalar field as color
            surf(obj.mesh.X, obj.mesh.Y, obj.mesh.Z, field);
            axis equal;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title(title_str);
            colormap('jet');
            colorbar;
            shading interp;
        end
    end
end
