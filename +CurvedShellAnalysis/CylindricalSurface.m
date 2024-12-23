classdef CylindricalSurface < CurvedShellAnalysis.Surface
    % CYLINDRICALSURFACE Class for cylindrical shell analysis
    
    properties
        R  % Radius of curvature (m)
        Angle = pi/2  % Angular span (radians)
        Direction = 'x'  % Cylinder axis direction ('x' or 'y')
    end
    
    methods
        function obj = CylindricalSurface(params)
            % Constructor for CylindricalSurface
            % params: struct with fields from Surface plus:
            %   R: radius
            %   Angle: angular span (optional, default: pi/2)
            %   Direction: axis direction (optional, default: 'x')
            
            obj = obj@CurvedShellAnalysis.Surface(params);
            obj.R = params.R;
            if isfield(params, 'Angle')
                obj.Angle = params.Angle;
            end
            if isfield(params, 'Direction')
                obj.Direction = params.Direction;
            end
            obj.validateParams();
        end
        
        function [X, Y, Z] = generateMesh(obj)
            % Generate cylindrical surface mesh
            if strcmpi(obj.Direction, 'x')
                [X, theta] = meshgrid(linspace(-obj.L/2, obj.L/2, obj.nx), ...
                                    linspace(-obj.Angle/2, obj.Angle/2, obj.ny));
                Y = obj.R * sin(theta);
                Z = obj.R * (1 - cos(theta));
            else  % y-direction
                [theta, Y] = meshgrid(linspace(-obj.Angle/2, obj.Angle/2, obj.nx), ...
                                    linspace(-obj.W/2, obj.W/2, obj.ny));
                X = obj.R * sin(theta);
                Z = obj.R * (1 - cos(theta));
            end
        end
        
        function validateParams(obj)
            % Validate parameters specific to cylindrical surface
            validateParams@CurvedShellAnalysis.Surface(obj);
            assert(obj.R > 0, 'Radius must be positive');
            assert(obj.Angle > 0 && obj.Angle <= 2*pi, ...
                  'Angle must be between 0 and 2π');
            assert(ismember(lower(obj.Direction), {'x','y'}), ...
                  'Direction must be ''x'' or ''y''');
        end
    end
end
