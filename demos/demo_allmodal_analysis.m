% 演示如何使用曲面壳模态分析包

% 清理工作空间
clear; clc; close all;

%% 1. 球面壳分析
fprintf('\n=== 球面壳分析 ===\n');

% 定义参数
params = struct();
params.L = 0.4;     % 长度 (m)
params.W = 0.3;     % 宽度 (m)
params.t = 0.002;   % 厚度 (m)
params.R = 0.5;     % 曲率半径 (m)
params.E = 2.1e11;  % 弹性模量 (Pa)
params.rho = 7800;  % 密度 (kg/m³)
params.nu = 0.3;    % 泊松比
params.zeta = 0.02; % 阻尼比
params.nx = 20;     % x方向节点数
params.ny = 15;     % y方向节点数

% 创建球面壳对象
spherical = CurvedShellAnalysis.SphericalSurface(params);

% 创建分析对象并运行分析
analysis = CurvedShellAnalysis.ModalAnalysis(spherical, 6);
analysis.analyze();

% 创建动态分析对象
dynamic = CurvedShellAnalysis.DynamicAnalysis(analysis);

% 设置中心点谐振力激励
center_node = ceil(params.nx * params.ny / 2);
dynamic.setHarmonicForce(center_node, 3, 1.0, 500);  % 1N, 500Hz

% 计算频率响应
freq_range = logspace(1, 3, 1000);  % 10Hz to 1kHz
harmonic_response = dynamic.harmonicResponse(freq_range);

% 计算瞬态响应
tspan = 0:0.0001:0.1;  % 0 to 0.1s
transient_response = dynamic.transientResponse(tspan);

% 创建热分析对象
thermal_params = struct();
thermal_params.alpha = 12e-6;  % 热膨胀系数 (1/K)
thermal_params.k = 45;         % 热导率 (W/m·K)
thermal_params.cp = 460;       % 比热容 (J/kg·K)

thermal = CurvedShellAnalysis.ThermalStructuralAnalysis(analysis, thermal_params);

% 设置温度场（示例：高斯分布）
temp_func = @(X,Y,Z) 25 + 75*exp(-(X.^2 + Y.^2)/(0.1^2));
thermal.setTemperatureField(temp_func);

% 计算热应力
[stress, strain] = thermal.calculateThermalStress();

% 创建屈曲分析对象
buckling = CurvedShellAnalysis.BucklingAnalysis(analysis);

% 分析不同载荷工况
% 1. 轴向压缩
buckling.setLoading('compression', 1e6);  % 1 MPa压缩
[lambda_comp, modes_comp] = buckling.linearBuckling(3);

% 2. 剪切
buckling.setLoading('shear', 5e5);  % 0.5 MPa剪切
[lambda_shear, modes_shear] = buckling.linearBuckling(3);

% 3. 外压
buckling.setLoading('pressure', 2e5);  % 0.2 MPa外压
[lambda_press, modes_press] = buckling.linearBuckling(3);

%% 2. 圆柱面壳分析
fprintf('\n=== 圆柱面壳分析 ===\n');

% 修改参数
params.R = 0.3;  % 圆柱半径
params.Angle = pi/2;  % 角度跨度
params.Direction = 'x';  % 轴向

% 创建圆柱面壳对象
cylindrical = CurvedShellAnalysis.CylindricalSurface(params);

% 创建新的分析对象并运行分析
analysis = CurvedShellAnalysis.ModalAnalysis(cylindrical, 6);
analysis.OutputFolder = 'ModalResults_Cylindrical';
analysis.analyze();

% 圆柱壳屈曲分析
buckling = CurvedShellAnalysis.BucklingAnalysis(analysis);
buckling.setLoading('compression', 1e6);
[lambda_cyl, modes_cyl] = buckling.linearBuckling(3);

%% 3. 圆锥面壳分析
fprintf('\n=== 圆锥面壳分析 ===\n');

% 修改参数
params = rmfield(params, {'R', 'Angle', 'Direction'});
params.R1 = 0.3;  % 底部半径
params.R2 = 0.1;  % 顶部半径

% 创建圆锥面壳对象
conical = CurvedShellAnalysis.ConicalSurface(params);

% 创建新的分析对象并运行分析
analysis = CurvedShellAnalysis.ModalAnalysis(conical, 6);
analysis.OutputFolder = 'ModalResults_Conical';
analysis.analyze();

% 圆锥壳热分析
thermal = CurvedShellAnalysis.ThermalStructuralAnalysis(analysis, thermal_params);
temp_func = @(X,Y,Z) 25 + 50*Z/max(Z(:));  % 线性温度梯度
thermal.setTemperatureField(temp_func);
[stress_cone, strain_cone] = thermal.calculateThermalStress();

%% 4. 椭球面壳分析
fprintf('\n=== 椭球面壳分析 ===\n');

% 修改参数
params = rmfield(params, {'R1', 'R2'});
params.Rx = 0.6;  % x方向曲率半径
params.Ry = 0.5;  % y方向曲率半径
params.Rz = 0.4;  % z方向曲率半径

% 创建椭球面壳对象
ellipsoidal = CurvedShellAnalysis.EllipsoidalSurface(params);

% 创建新的分析对象并运行分析
analysis = CurvedShellAnalysis.ModalAnalysis(ellipsoidal, 6);
analysis.OutputFolder = 'ModalResults_Ellipsoidal';
analysis.analyze();

% 椭球壳屈曲分析
buckling = CurvedShellAnalysis.BucklingAnalysis(analysis);
buckling.setLoading('pressure', 2e5);
[lambda_ell, modes_ell] = buckling.linearBuckling(3);

%% 5. 自定义曲面分析
fprintf('\n=== 自定义曲面分析 ===\n');

% 定义自定义曲面函数（例如：双曲抛物面）
params.SurfaceFunction = @(X,Y) 0.1 * (X.^2/params.Rx^2 - Y.^2/params.Ry^2);

% 创建自定义曲面对象
custom = CurvedShellAnalysis.CustomSurface(params);

% 创建新的分析对象并运行分析
analysis = CurvedShellAnalysis.ModalAnalysis(custom, 6);
analysis.OutputFolder = 'ModalResults_Custom';
analysis.analyze();

% 自定义曲面热-结构耦合分析
thermal = CurvedShellAnalysis.ThermalStructuralAnalysis(analysis, thermal_params);
temp_func = @(X,Y,Z) 25 + 25*sin(2*pi*X/params.L).*cos(2*pi*Y/params.W);
thermal.setTemperatureField(temp_func);
[stress_custom, strain_custom] = thermal.calculateThermalStress();

fprintf('\n所有分析完成！\n');

% 显示各种曲面的临界屈曲载荷比较
fprintf('\n临界屈曲载荷比较:\n');
fprintf('球壳 - 压缩: %.2e Pa\n', min(lambda_comp));
fprintf('球壳 - 剪切: %.2e Pa\n', min(lambda_shear));
fprintf('球壳 - 外压: %.2e Pa\n', min(lambda_press));
fprintf('圆柱壳 - 压缩: %.2e Pa\n', min(lambda_cyl));
fprintf('椭球壳 - 外压: %.2e Pa\n', min(lambda_ell));
