classdef CustomSurface < CurvedShellAnalysis.Surface
    % CUSTOMSURFACE Class for user-defined curved surface analysis
    
    properties
        SurfaceFunction  % Function handle for Z = f(X,Y)
    end
    
    methods
        function obj = CustomSurface(params)
            % Constructor for CustomSurface
            % params: struct with fields from Surface plus SurfaceFunction
            
            obj = obj@CurvedShellAnalysis.Surface(params);
            obj.SurfaceFunction = params.SurfaceFunction;
            obj.validateParams();
        end
        
        function [X, Y, Z] = generateMesh(obj)
            % Generate custom surface mesh using the provided function
            [X, Y] = meshgrid(linspace(-obj.L/2, obj.L/2, obj.nx), ...
                            linspace(-obj.W/2, obj.W/2, obj.ny));
            
            % Calculate Z coordinates using the custom function
            Z = obj.SurfaceFunction(X, Y);
        end
        
        function validateParams(obj)
            % Validate parameters specific to custom surface
            validateParams@CurvedShellAnalysis.Surface(obj);
            assert(isa(obj.SurfaceFunction, 'function_handle'), ...
                  'SurfaceFunction must be a function handle');
            
            % Test the function with a simple input
            try
                [X, Y] = meshgrid([-1 1], [-1 1]);
                Z = obj.SurfaceFunction(X, Y);
                assert(isnumeric(Z) && size(Z,1) == 2 && size(Z,2) == 2, ...
                      'SurfaceFunction must return a numeric array of same size as input');
            catch
                error('Invalid SurfaceFunction: must accept X,Y matrices and return Z matrix');
            end
        end
    end
end
