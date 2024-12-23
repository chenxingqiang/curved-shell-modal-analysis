classdef Surface
    % SURFACE Base class for curved surfaces
    
    properties
        L  % Length in x direction (m)
        W  % Width in y direction (m)
        t  % Thickness (m)
        E  % Young's modulus (Pa)
        rho  % Density (kg/m³)
        nu  % Poisson's ratio
        zeta  % Damping ratio
        nx  % Number of nodes in x direction
        ny  % Number of nodes in y direction
    end
    
    methods
        function obj = Surface(params)
            % Constructor for Surface class
            % params: struct with fields L, W, t, E, rho, nu, zeta, nx, ny
            
            if nargin > 0
                obj.L = params.L;
                obj.W = params.W;
                obj.t = params.t;
                obj.E = params.E;
                obj.rho = params.rho;
                obj.nu = params.nu;
                obj.zeta = params.zeta;
                obj.nx = params.nx;
                obj.ny = params.ny;
            end
        end
        
        function [X, Y, Z] = generateMesh(obj)
            % Generate basic mesh (to be overridden by subclasses)
            [X, Y] = meshgrid(linspace(-obj.L/2, obj.L/2, obj.nx), ...
                            linspace(-obj.W/2, obj.W/2, obj.ny));
            Z = zeros(size(X));  % Flat surface by default
        end
        
        function validateParams(obj)
            % Validate parameters
            assert(obj.L > 0, 'Length must be positive');
            assert(obj.W > 0, 'Width must be positive');
            assert(obj.t > 0, 'Thickness must be positive');
            assert(obj.E > 0, 'Young''s modulus must be positive');
            assert(obj.rho > 0, 'Density must be positive');
            assert(obj.nu >= 0 && obj.nu < 0.5, 'Poisson''s ratio must be between 0 and 0.5');
            assert(obj.zeta >= 0, 'Damping ratio must be non-negative');
            assert(obj.nx > 1, 'Number of nodes in x direction must be greater than 1');
            assert(obj.ny > 1, 'Number of nodes in y direction must be greater than 1');
        end
    end
end
