classdef ConicalSurface < CurvedShellAnalysis.Surface
    % CONICALSURFACE Class for conical shell analysis
    
    properties
        Alpha  % Cone semi-angle (radians)
        R1     % Bottom radius (m)
        R2     % Top radius (m)
    end
    
    methods
        function obj = ConicalSurface(params)
            % Constructor for ConicalSurface
            % params: struct with fields from Surface plus:
            %   Alpha: cone semi-angle or
            %   R1, R2: bottom and top radii
            
            obj = obj@CurvedShellAnalysis.Surface(params);
            if isfield(params, 'Alpha')
                obj.Alpha = params.Alpha;
                obj.R1 = params.R1;
                obj.R2 = obj.R1 - obj.L*tan(obj.Alpha);
            else
                obj.R1 = params.R1;
                obj.R2 = params.R2;
                obj.Alpha = atan2(obj.R1 - obj.R2, obj.L);
            end
            obj.validateParams();
            obj.createMesh();
        end
        
        function createMesh(obj)
            % Create conical surface mesh
            [theta, h] = meshgrid(linspace(0, 2*pi, obj.nx), ...
                                linspace(0, obj.L, obj.ny));
            
            % Radius varies linearly with height
            r = obj.R1 - (obj.R1 - obj.R2) * h/obj.L;
            
            % Store mesh data
            obj.mesh = struct();
            obj.mesh.X = r .* cos(theta);
            obj.mesh.Y = r .* sin(theta);
            obj.mesh.Z = h;
        end
        
        function validateParams(obj)
            % Validate parameters specific to conical surface
            if ~isempty(obj.R1) && ~isempty(obj.R2)
                assert(obj.R1 > 0 && obj.R2 >= 0, 'Radii must be non-negative');
                assert(obj.R1 >= obj.R2, 'Bottom radius must be larger than top radius');
            end
            if ~isempty(obj.Alpha)
                assert(obj.Alpha > 0 && obj.Alpha < pi/2, ...
                      'Cone angle must be between 0 and π/2');
            end
            validateParams@CurvedShellAnalysis.Surface(obj);
        end
    end
end
