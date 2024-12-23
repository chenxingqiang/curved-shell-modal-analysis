classdef DynamicAnalysis < handle
    % DYNAMICANALYSIS Class for dynamic analysis of curved surfaces
    
    properties
        Surface  % Surface object
        ModalAnalysis  % Reference to modal analysis object
        TimeSpan  % Time vector for transient analysis
        Forces  % Structure containing force information
        OutputFolder = 'DynamicResults'  % Output folder
    end
    
    methods
        function obj = DynamicAnalysis(modalAnalysis)
            % Constructor
            obj.ModalAnalysis = modalAnalysis;
            obj.Surface = modalAnalysis.Surface;
            
            % Create output folder
            if ~exist(obj.OutputFolder, 'dir')
                mkdir(obj.OutputFolder);
            end
        end
        
        function setHarmonicForce(obj, node, direction, amplitude, frequency)
            % Set harmonic force parameters
            obj.Forces.Type = 'harmonic';
            obj.Forces.Node = node;
            obj.Forces.Direction = direction;
            obj.Forces.Amplitude = amplitude;
            obj.Forces.Frequency = frequency;
        end
        
        function setImpulseForce(obj, node, direction, amplitude)
            % Set impulse force parameters
            obj.Forces.Type = 'impulse';
            obj.Forces.Node = node;
            obj.Forces.Direction = direction;
            obj.Forces.Amplitude = amplitude;
        end
        
        function setRandomForce(obj, node, direction, PSD)
            % Set random force parameters
            obj.Forces.Type = 'random';
            obj.Forces.Node = node;
            obj.Forces.Direction = direction;
            obj.Forces.PSD = PSD;  % Power Spectral Density function
        end
        
        function response = harmonicResponse(obj, freq_range)
            % Calculate harmonic response
            fprintf('计算谐响应...\n');
            
            w = 2*pi*freq_range;
            ndof = obj.ModalAnalysis.DOFsPerNode;
            force_dof = (obj.Forces.Node-1)*ndof + obj.Forces.Direction;
            
            % Initialize response matrix
            response.Frequency = freq_range;
            response.Amplitude = zeros(size(freq_range));
            response.Phase = zeros(size(freq_range));
            
            % Calculate response at each frequency
            for i = 1:length(w)
                % Dynamic stiffness matrix
                Z = -w(i)^2*obj.ModalAnalysis.MassMatrix + ...
                    1i*w(i)*obj.ModalAnalysis.DampingMatrix + ...
                    obj.ModalAnalysis.StiffnessMatrix;
                
                % Force vector
                F = zeros(size(Z,1), 1);
                F(force_dof) = obj.Forces.Amplitude;
                
                % Solve for displacement
                U = Z\F;
                
                % Store amplitude and phase
                response.Amplitude(i) = abs(U(force_dof));
                response.Phase(i) = angle(U(force_dof));
            end
            
            % Plot results
            obj.plotHarmonicResponse(response);
        end
        
        function response = transientResponse(obj, tspan)
            % Calculate transient response
            fprintf('计算瞬态响应...\n');
            
            obj.TimeSpan = tspan;
            ndof = obj.ModalAnalysis.DOFsPerNode;
            force_dof = (obj.Forces.Node-1)*ndof + obj.Forces.Direction;
            
            % Mass, damping, and stiffness matrices
            M = obj.ModalAnalysis.MassMatrix;
            C = obj.ModalAnalysis.DampingMatrix;
            K = obj.ModalAnalysis.StiffnessMatrix;
            
            % Convert to state-space form
            n = size(M,1);
            A = [zeros(n), eye(n); 
                 -M\K, -M\C];
            B = [zeros(n,1); M\[zeros(force_dof-1,1); 1; zeros(n-force_dof,1)]];
            
            % Initial conditions
            x0 = zeros(2*n,1);
            
            % Time integration using ode45
            [t, x] = ode45(@(t,x) obj.stateSpaceModel(t,x,A,B), tspan, x0);
            
            % Extract displacement response
            response.Time = t;
            response.Displacement = x(:,force_dof);
            response.Velocity = x(:,n+force_dof);
            
            % Plot results
            obj.plotTransientResponse(response);
        end
        
        function response = randomResponse(obj, freq_range)
            % Calculate random vibration response
            fprintf('计算随机响应...\n');
            
            w = 2*pi*freq_range;
            ndof = obj.ModalAnalysis.DOFsPerNode;
            force_dof = (obj.Forces.Node-1)*ndof + obj.Forces.Direction;
            
            % Initialize PSD matrix
            response.Frequency = freq_range;
            response.PSD = zeros(size(freq_range));
            
            % Calculate PSD at each frequency
            for i = 1:length(w)
                % Frequency response function
                Z = -w(i)^2*obj.ModalAnalysis.MassMatrix + ...
                    1i*w(i)*obj.ModalAnalysis.DampingMatrix + ...
                    obj.ModalAnalysis.StiffnessMatrix;
                H = inv(Z);
                
                % Calculate response PSD
                response.PSD(i) = abs(H(force_dof,force_dof))^2 * ...
                                obj.Forces.PSD(freq_range(i));
            end
            
            % Calculate RMS response
            df = freq_range(2) - freq_range(1);
            response.RMS = sqrt(sum(response.PSD)*df);
            
            % Plot results
            obj.plotRandomResponse(response);
        end
        
        function dxdt = stateSpaceModel(obj, t, x, A, B)
            % State-space model for time integration
            if strcmpi(obj.Forces.Type, 'harmonic')
                f = obj.Forces.Amplitude * ...
                    sin(2*pi*obj.Forces.Frequency*t);
            elseif strcmpi(obj.Forces.Type, 'impulse')
                f = obj.Forces.Amplitude * ...
                    (abs(t) < 1e-3);  % Approximate impulse
            else
                f = 0;
            end
            dxdt = A*x + B*f;
        end
        
        function plotHarmonicResponse(obj, response)
            % Plot harmonic response
            figure('Position', [100 100 800 600]);
            
            % Magnitude plot
            subplot(2,1,1);
            semilogx(response.Frequency, 20*log10(response.Amplitude), ...
                     'LineWidth', 2);
            grid on;
            title('频率响应幅值');
            xlabel('频率 (Hz)');
            ylabel('幅值 (dB)');
            
            % Phase plot
            subplot(2,1,2);
            semilogx(response.Frequency, response.Phase*180/pi, ...
                     'LineWidth', 2);
            grid on;
            title('频率响应相位');
            xlabel('频率 (Hz)');
            ylabel('相位 (度)');
            
            % Save figure
            saveas(gcf, fullfile(obj.OutputFolder, 'HarmonicResponse.png'));
        end
        
        function plotTransientResponse(obj, response)
            % Plot transient response
            figure('Position', [100 100 800 600]);
            
            % Displacement plot
            subplot(2,1,1);
            plot(response.Time, response.Displacement, 'LineWidth', 2);
            grid on;
            title('位移响应');
            xlabel('时间 (s)');
            ylabel('位移 (m)');
            
            % Velocity plot
            subplot(2,1,2);
            plot(response.Time, response.Velocity, 'LineWidth', 2);
            grid on;
            title('速度响应');
            xlabel('时间 (s)');
            ylabel('速度 (m/s)');
            
            % Save figure
            saveas(gcf, fullfile(obj.OutputFolder, 'TransientResponse.png'));
        end
        
        function plotRandomResponse(obj, response)
            % Plot random response
            figure('Position', [100 100 800 400]);
            
            loglog(response.Frequency, response.PSD, 'LineWidth', 2);
            grid on;
            title(sprintf('随机响应功率谱密度 (RMS = %.3e)', response.RMS));
            xlabel('频率 (Hz)');
            ylabel('PSD (m²/Hz)');
            
            % Save figure
            saveas(gcf, fullfile(obj.OutputFolder, 'RandomResponse.png'));
        end
    end
end
