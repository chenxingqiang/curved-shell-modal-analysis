%% 曲面板模态分析 - 增强版本
% Author: Xingqiang Chen
% Date: 2024-12-23

% 清理工作空间
clear; clc; close all;

% 创建结果文件夹
if ~exist('results', 'dir')
    mkdir('results');
end
if ~exist('../docs/images', 'dir')
    mkdir('../docs/images');
end

%% 定义几何参数
% 长度 (m)
L = 0.4;    
% 宽度 (m)
W = 0.3;    
% 厚度 (m)
t = 0.002;  
% 曲率半径 (m)
R = 0.5;    

%% 定义材料属性（钢）
% 弹性模量 (Pa)
E = 2.1e11;  
% 密度 (kg/m³)
rho = 7800;  
% 泊松比
nu = 0.3;    
% 阻尼比
zeta = 0.02; 

%% 创建网格
% x方向的节点数
nx = 20;  
% y方向的节点数
ny = 15;  
[X, Y] = meshgrid(linspace(-L/2, L/2, nx), linspace(-W/2, W/2, ny));
% 计算曲面高度
Z = R - sqrt(R^2 - Y.^2);  

%% 创建有限元模型
% 每个节点的自由度数：u, v, w, θx, θy
ndof = 5;  
% 总自由度数
N = nx * ny * ndof;  
% 质量矩阵
M = sparse(N, N);  
% 刚度矩阵
K = sparse(N, N);  
% 阻尼矩阵
C = sparse(N, N);  

%% 定义单元参数
% x方向单元长度
dx = L/(nx-1);
% y方向单元长度
dy = W/(ny-1);
% 弯曲刚度
D = E*t^3/(12*(1-nu^2));  
% 剪切模量
G = E/(2*(1+nu));         
% 剪切修正系数
kappa = 5/6;              

%% 创建进度条
fprintf('正在组装全局矩阵...\n');
h = waitbar(0, '组装单元矩阵...');

%% 组装全局矩阵
for i = 1:nx-1
    for j = 1:ny-1
        % 更新进度条
        waitbar(((i-1)*(ny-1) + j)/((nx-1)*(ny-1)), h);
        
        % 单元节点编号
        n1 = j + (i-1)*ny;
        n2 = j + i*ny;
        n3 = j+1 + i*ny;
        n4 = j+1 + (i-1)*ny;
        
        % 获取节点坐标
        x1 = X(j,i); y1 = Y(j,i); z1 = Z(j,i);
        x2 = X(j,i+1); y2 = Y(j,i+1); z2 = Z(j,i+1);
        x3 = X(j+1,i+1); y3 = Y(j+1,i+1); z3 = Z(j+1,i+1);
        x4 = X(j+1,i); y4 = Y(j+1,i); z4 = Z(j+1,i);
        
        % 计算单元法向量
        v1 = [x2-x1; y2-y1; z2-z1];
        v2 = [x4-x1; y4-y1; z4-z1];
        n = cross(v1, v2);
        n = n/norm(n);
        
        % 单元自由度
        dofs = [];
        for n = [n1 n2 n3 n4]
            dofs = [dofs, (n-1)*ndof+1:(n-1)*ndof+5];
        end
        
        % 计算单元矩阵
        A = dx*dy;  % 单元面积
        
        % 简化的单元刚度矩阵
        ke = zeros(20, 20);
        
        % 膜变形部分
        ke(1:12, 1:12) = E*t/(1-nu^2) * A/4 * [
            4  2  0  2  1  0  1  0  0  2  1  0
            2  4  0  1  2  0  0  1  0  1  2  0
            0  0  1  0  0  0  0  0  0  0  0  0
            2  1  0  4  2  0  2  1  0  1  0  0
            1  2  0  2  4  0  1  2  0  0  1  0
            0  0  0  0  0  1  0  0  0  0  0  0
            1  0  0  2  1  0  4  2  0  2  1  0
            0  1  0  1  2  0  2  4  0  1  2  0
            0  0  0  0  0  0  0  0  1  0  0  0
            2  1  0  1  0  0  2  1  0  4  2  0
            1  2  0  0  1  0  1  2  0  2  4  0
            0  0  0  0  0  0  0  0  0  0  0  1
        ];
        
        % 弯曲部分
        ke(13:20, 13:20) = D*A/4 * [
            4  1  1  0  1  0  0  0
            1  4  0  1  0  1  0  0
            1  0  4  1  0  0  1  0
            0  1  1  4  0  0  0  1
            1  0  0  0  4  1  1  0
            0  1  0  0  1  4  0  1
            0  0  1  0  1  0  4  1
            0  0  0  1  0  1  1  4
        ];
        
        % 质量矩阵
        me = zeros(20, 20);
        
        % 平动质量
        me(1:12, 1:12) = rho*t*A/36 * [
            4  0  0  2  0  0  1  0  0  2  0  0
            0  4  0  0  2  0  0  1  0  0  2  0
            0  0  4  0  0  2  0  0  1  0  0  2
            2  0  0  4  0  0  2  0  0  1  0  0
            0  2  0  0  4  0  0  2  0  0  1  0
            0  0  2  0  0  4  0  0  2  0  0  1
            1  0  0  2  0  0  4  0  0  2  0  0
            0  1  0  0  2  0  0  4  0  0  2  0
            0  0  1  0  0  2  0  0  4  0  0  2
            2  0  0  1  0  0  2  0  0  4  0  0
            0  2  0  0  1  0  0  2  0  0  4  0
            0  0  2  0  0  1  0  0  2  0  0  4
        ];
        
        % 转动质量
        me(13:20, 13:20) = rho*t^3*A/420 * [
            16  0 -8   0  4  0 -8   0
             0 16  0  -8  0  4  0  -8
            -8  0 16   0 -8  0  4   0
             0 -8  0  16  0 -8  0   4
             4  0 -8   0 16  0 -8   0
             0  4  0  -8  0 16  0  -8
            -8  0  4   0 -8  0 16   0
             0 -8  0   4  0 -8  0  16
        ];
        
        % 组装到全局矩阵
        K(dofs,dofs) = K(dofs,dofs) + ke;
        M(dofs,dofs) = M(dofs,dofs) + me;
    end
end

%% 关闭进度条
close(h);

%% 应用边界条件（固定两端）
fprintf('应用边界条件...\n');
fixed_dofs = [];
for i = [1, nx]
    for j = 1:ny
        n = j + (i-1)*ny;
        fixed_dofs = [fixed_dofs, (n-1)*ndof+1:(n-1)*ndof+5];
    end
end
free_dofs = setdiff(1:N, fixed_dofs);

%% 求解特征值问题
fprintf('求解特征值问题...\n');
[V, D] = eig(full(K(free_dofs,free_dofs)), full(M(free_dofs,free_dofs)));
[frequencies, idx] = sort(sqrt(real(diag(D)))/(2*pi));  % 确保频率为实数并按升序排列
V = V(:,idx);  % 对应排序特征向量
frequencies = frequencies(1:6);  % 只取前6阶
V = V(:,1:6);  % 只取前6阶模态

%% 构建完整模态
modes = zeros(N, 6);
for i = 1:6
    mode_free = V(:,i);
    mode_full = zeros(N, 1);
    mode_full(free_dofs) = mode_free;
    modes(:,i) = mode_full;
end

%% 计算阻尼矩阵（Rayleigh阻尼）
alpha = 2*zeta*sqrt(frequencies(1)*frequencies(6));  % 质量比例系数
beta = 2*zeta/(sqrt(frequencies(1)*frequencies(6))); % 刚度比例系数
C = alpha*M + beta*K;

%% 保存频率结果
fid = fopen(fullfile('results', 'NaturalFrequencies.txt'), 'w');
fprintf(fid, '曲面板模态分析结果\n\n');
fprintf(fid, '几何参数:\n');
fprintf(fid, '长度: %.3f m\n', L);
fprintf(fid, '宽度: %.3f m\n', W);
fprintf(fid, '厚度: %.3f m\n', t);
fprintf(fid, '曲率半径: %.3f m\n\n', R);
fprintf(fid, '材料属性:\n');
fprintf(fid, '弹性模量: %.2e Pa\n', E);
fprintf(fid, '密度: %.1f kg/m³\n', rho);
fprintf(fid, '泊松比: %.2f\n', nu);
fprintf(fid, '阻尼比: %.3f\n\n', zeta);
fprintf(fid, '前6阶固有频率:\n');
for i = 1:length(frequencies)
    fprintf(fid, '第%d阶: %.2f Hz\n', i, frequencies(i));
end
fprintf(fid, '\n阻尼参数:\n');
fprintf(fid, 'α (质量比例系数): %.6e\n', alpha);
fprintf(fid, 'β (刚度比例系数): %.6e\n', beta);
fclose(fid);

%% 绘制模态振型
fprintf('绘制模态振型...\n');
figure('Name', 'Modal Analysis', 'Position', [100 100 1200 800]);
for i = 1:6
    subplot(2,3,i);
    % 提取位移分量
    w = zeros(ny, nx);
    for j = 1:nx
        for k = 1:ny
            n = k + (j-1)*ny;
            w(k,j) = modes((n-1)*ndof+3, i);  % 取出z方向位移
        end
    end
    
    % 归一化位移
    w = w/max(abs(w(:)));
    
    % 创建变形后的网格
    Zdef = Z + w*0.1*max(abs(Z(:)));  % 放大变形以便观察
    
    % 绘制3D曲面
    surf(X, Y, Zdef);
    colormap('jet');
    shading interp;
    colorbar;
    title(sprintf('Mode %d: %.2f Hz', i, frequencies(i)));
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    axis equal;
    view(45, 30);
end
saveas(gcf, fullfile('results', 'ModalAnalysis.png'));
saveas(gcf, fullfile('../docs/images', 'ModalAnalysis.png'));

%% 绘制频率响应函数
fprintf('计算频率响应函数...\n');
f = logspace(0, 4, 1000);  % 0.1 Hz到10 kHz
w = 2*pi*f;
H = zeros(length(f), 1);

% 选择观察点（面板中心）
obs_node = ceil(nx*ny/2);
obs_dof = (obs_node-1)*ndof + 3;  % z方向自由度

% 计算FRF
for i = 1:length(f)
    Z = -w(i)^2*M + 1i*w(i)*C + K;
    Z = Z(free_dofs,free_dofs);
    F = zeros(length(free_dofs), 1);
    F(obs_dof) = 1;  % 单位力激励
    U = Z\F;
    H(i) = abs(U(obs_dof));
end

% 绘制FRF
figure('Name', 'Frequency Response', 'Position', [100 100 800 600]);
semilogx(f, 20*log10(H), 'LineWidth', 2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Frequency Response Function');
saveas(gcf, fullfile('results', 'FrequencyResponse.png'));
saveas(gcf, fullfile('../docs/images', 'FrequencyResponse.png'));

%% 保存所有图像
figHandles = findall(0, 'Type', 'figure');
for i = 1:length(figHandles)
    fig = figHandles(i);
    figName = get(fig, 'Name');
    if ~isempty(figName)
        % 转换空格为下划线并删除特殊字符
        figName = lower(regexprep(figName, '\s+', '_'));
        figName = regexprep(figName, '[^a-z0-9_]', '');
        
        % 保存到results和docs/images
        print(fig, fullfile('results', figName), '-dpng', '-r300');
        print(fig, fullfile('../docs/images', figName), '-dpng', '-r300');
    end
end

fprintf('\n分析完成！结果已保存到 results/ 和 docs/images/ 目录\n');