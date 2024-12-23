classdef CylindricalSurface < CurvedShellAnalysis.Surface
    % CylindricalSurface Class for cylindrical shell analysis
    
    properties
        R       % Radius
        zeta    % Damping ratio
    end
    
    methods
        function obj = CylindricalSurface(params)
            % Constructor
            obj = obj@CurvedShellAnalysis.Surface(params);
            
            % Set cylindrical specific parameters
            if isfield(params, 'R')
                obj.R = params.R;
            else
                error('Radius (R) must be specified for cylindrical surface');
            end
            
            if isfield(params, 'zeta')
                obj.zeta = params.zeta;
            else
                obj.zeta = 0.02;  % Default damping ratio
            end
            
            % Validate and create mesh
            obj.validateParams();
            obj.createMesh();
        end
        
        function createMesh(obj)
            % Create mesh for cylindrical surface
            [x, y] = meshgrid(linspace(0, obj.L, obj.nx), ...
                            linspace(-obj.W/2, obj.W/2, obj.ny));
            
            % Calculate z coordinates for cylindrical surface
            theta = y / obj.R;  % Angular position
            z = obj.R * (1 - cos(theta));
            
            % Store mesh data
            obj.mesh = struct('X', x, 'Y', y, 'Z', z);
        end
        
        function validateParams(obj)
            % Validate parameters specific to cylindrical surface
            if ~isempty(obj.R)
                assert(obj.R > 0, 'Radius must be positive');
            end
            if ~isempty(obj.zeta)
                assert(obj.zeta >= 0 && obj.zeta < 1, ...
                      'Damping ratio must be between 0 and 1');
            end
            validateParams@CurvedShellAnalysis.Surface(obj);
        end
    end
end
