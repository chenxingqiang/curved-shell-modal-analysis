classdef SphericalSurface < CurvedShellAnalysis.Surface
    % SPHERICALSURFACE Class for spherical shell analysis
    
    properties
        R  % Radius of curvature (m)
    end
    
    methods
        function obj = SphericalSurface(params)
            % Constructor for SphericalSurface
            % params: struct with fields from Surface plus R (radius)
            
            obj = obj@CurvedShellAnalysis.Surface(params);
            obj.R = params.R;
            obj.validateParams();
        end
        
        function [X, Y, Z] = generateMesh(obj)
            % Generate spherical surface mesh
            [X, Y] = meshgrid(linspace(-obj.L/2, obj.L/2, obj.nx), ...
                            linspace(-obj.W/2, obj.W/2, obj.ny));
            
            % Calculate Z coordinates for spherical surface
            Z = obj.R - sqrt(obj.R^2 - (X.^2 + Y.^2));
        end
        
        function validateParams(obj)
            % Validate parameters specific to spherical surface
            validateParams@CurvedShellAnalysis.Surface(obj);
            assert(obj.R > max(obj.L/2, obj.W/2), ...
                  'Radius must be larger than half the maximum dimension');
        end
    end
end
