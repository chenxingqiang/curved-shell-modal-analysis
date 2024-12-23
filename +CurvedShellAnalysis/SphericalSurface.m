classdef SphericalSurface < CurvedShellAnalysis.Surface
    % SphericalSurface Class for spherical shell surfaces
    
    properties
        R           % Radius
        ThetaSpan   % Polar angle span
        PhiSpan     % Azimuthal angle span
    end
    
    methods
        function obj = SphericalSurface(params)
            % Constructor for SphericalSurface
            obj = obj@CurvedShellAnalysis.Surface(params);
            
            % Set spherical specific parameters
            if isfield(params, 'R')
                obj.R = params.R;
            else
                error('Radius (R) must be specified for spherical surface');
            end
            
            if isfield(params, 'phi')
                obj.PhiSpan = params.phi;
            else
                obj.PhiSpan = pi/2;  % Default to quarter sphere
            end
            
            if isfield(params, 'theta')
                obj.ThetaSpan = params.theta;
            else
                obj.ThetaSpan = pi/2;  % Default to quarter sphere
            end
            
            % Now validate all parameters
            obj.validateParams();
            
            % Create mesh
            obj.createMesh();
        end
        
        function createMesh(obj)
            % Create spherical mesh
            [phi, theta] = meshgrid(linspace(0, obj.PhiSpan, obj.nx), ...
                                  linspace(0, obj.ThetaSpan, obj.ny));
            
            obj.mesh = struct();
            obj.mesh.X = obj.R * sin(theta) .* cos(phi);
            obj.mesh.Y = obj.R * sin(theta) .* sin(phi);
            obj.mesh.Z = obj.R * cos(theta);
        end
        
        function plotSurface(obj)
            % Plot the spherical surface
            surf(obj.mesh.X, obj.mesh.Y, obj.mesh.Z);
            axis equal;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title('Spherical Surface');
            colormap('jet');
            colorbar;
        end
        
        function validateParams(obj)
            % Validate parameters specific to spherical surface
            if ~isempty(obj.R)
                assert(obj.R > 0, 'Radius must be positive');
            end
            if ~isempty(obj.PhiSpan)
                assert(obj.PhiSpan > 0 && obj.PhiSpan <= 2*pi, ...
                      'Phi span must be between 0 and 2π');
            end
            if ~isempty(obj.ThetaSpan)
                assert(obj.ThetaSpan > 0 && obj.ThetaSpan <= pi, ...
                      'Theta span must be between 0 and π');
            end
            validateParams@CurvedShellAnalysis.Surface(obj);
        end
    end
end
