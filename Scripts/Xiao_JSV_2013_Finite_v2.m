% 目的：找到FRF求解方法
% 文献：Flexural wave propagation in beams with periodically attached vibration absorbers: Band-gap behavior and band formation mechanisms
% 
% 状态：失败，不是波数k的问题。k与fre呈二次函数关系，正常。
% 备注：尝试更换k的计算方法。
% ================================================
%%
%[text] ## 参数设置
clear;clc


%% 1. 几何参数 (Table 1)
geometry.h = 0.015;      % 梁高度 [m]
geometry.l = 0.01;       % 间距 [m]
geometry.r1 = 0.013;     % 内半径 [m]
geometry.r2 = 0.015;     % 中间半径 [m]
geometry.r3 = 0.018;     % 外半径 [m]
geometry.L = 0.6;        % 梁总长度 [m]
geometry.N = 6;         % 单元个数

%% 2. 材料参数 (Table 2)
% 铝层参数
material.aluminum.E = 7e10;          % 杨氏模量 [Pa]
material.aluminum.G = 2.7e10;        % 剪切模量 [Pa]
material.aluminum.rho = 2700;        % 密度 [kg/m^3]

% 橡胶层参数
material.rubber.E = 7.748e5;         % 杨氏模量 [Pa]
material.rubber.G = 2.6e5;           % 剪切模量 [Pa]
material.rubber.rho = 1300;          % 密度 [kg/m^3]

% 铜层参数
material.copper.E = 1.65e11;         % 杨氏模量 [Pa]
material.copper.G = 7.53e10;         % 剪切模量 [Pa]
material.copper.rho = 8700;          % 密度 [kg/m^3]

%% 3. 截面特性计算
% 梁的横截面积
geometry.A = pi*(geometry.r3^2-geometry.r2^2);
% 梁的惯性矩
geometry.I = pi*(geometry.r3^4-geometry.r2^4)/4;

%% 4. DVA参数计算
% 形状参数 S = l / [(r1 + r2) * ln(r2/r1)]
geometry.DVA.S = geometry.l / ((geometry.r1 + geometry.r2) * log(geometry.r2/geometry.r1));

% DVA质量 m = ρ_Cu * π * （r3²-r2²） * h
material.DVA.m = material.copper.rho * pi * (geometry.r1^2) * geometry.h;

% DVA弹簧刚度 k = [π(5 + 3.29S²)Gl] / ln(r₂/r₁)
material.DVA.k = (pi * (5 + 3.29*geometry.DVA.S^2) * material.rubber.G * geometry.l) / ...
             log(geometry.r2/geometry.r1);






%%
%[text] ## 计算和结果

%% 频率扫描范围
omega_values = linspace(0.5, 12000, 12000);
fre_values = omega_values ./ (2*pi());  % 转换为Hz
T_N_results = zeros(size(omega_values)); % 存储透射率结果

% 检查波束随频率变化，后续画图观察
k_record = zeros(size(omega_values));

%% 主循环
for i = 1:length(omega_values)
    omega = omega_values(i);
    
    % 计算当前频率下的透射率
    [T_N, k] = calculate_transmission_rate(omega, material, geometry);
    % 转换为分贝(dB)：T_N(dB) = 20*log10(|T_N|)
    % T_N_results(i) = T_N;
    T_N_results(i) = 20 * log10(abs(T_N));
    
    k_record(i) = k;
end

%% 绘制结果
figure; %[output:599fbcb9]
plot(fre_values, T_N_results, 'b-', 'LineWidth', 2); %[output:599fbcb9]
xlabel('频率 [Hz]'); %[output:599fbcb9]
ylabel('透射率 T_N [dB]'); %[output:599fbcb9]
title('透射率随频率变化曲线'); %[output:599fbcb9]
grid on; %[output:599fbcb9]
xlim([0 2000]); %[output:599fbcb9]


plot_k(fre_values, k_record) %[output:6f16dfd2]
%%
function [T_N, k_record] = calculate_transmission_rate(omega, material, geometry)
    % 构建全局矩阵系统
    % 输入参数：
    %   material : 材料参数结构体
    %   geometry : 几何参数结构体
    %   omega : 角频率 [rad/s]
    %
    % 输出参数：
    %   T_N : 透射率




    % 计算DVA动态刚度
    D0 = calculate_DVA_stiffness(omega, material);

    % 单元长度
    geometry.ucL = geometry.L / geometry.N;
    % La和Lb段单元长度
    La = geometry.ucL/2;  % 假设La和Lb各占一半长度
    Lb = geometry.ucL/2;

    % 计算单元矩阵（La和Lb段）
    [D_La,k_record] = spectral_element_matrix(omega, material, geometry, La);
    D_Lb = spectral_element_matrix(omega, material, geometry, Lb);
    
    % 提取子矩阵块
    D_LL = D_La(1:2,1:2);       % D11(La)
    D_LI = D_La(1:2,3:4);       % D12(La)
    D_IL = D_La(3:4,1:2);       % D21(La)
    D_II_La = D_La(3:4,3:4);    % D22(La)
    
    D_II_Lb = D_Lb(1:2,1:2);    % D11(Lb)
    D_IR = D_Lb(1:2,3:4);       % D12(Lb)
    D_RI = D_Lb(3:4,1:2);       % D21(Lb)
    D_RR = D_Lb(3:4,3:4);       % D22(Lb)

    % 凝聚矩阵
    % 计算公共逆矩阵项
    D_II = D_II_La + D_II_Lb + [D0 0; 0 0];
    inv_D_II = inv(D_II);

    D_LL_ = D_LL - D_LI * inv_D_II * D_IL;
    D_LR_ = -D_LI * inv_D_II * D_IR;
    D_RL_ = -D_RI * inv_D_II * D_IL;
    D_RR_ = D_RR - D_RI * inv_D_II * D_IR;

    % 初始化全局矩阵 (n_elements+1节点，每个节点有2个自由度)
    global_dof = 2 * (geometry.N + 1);
    D_global = zeros(global_dof, global_dof);

    % 激励和响应位置设置
    excitation_dof = 1;  % 激励节点
    response_dof = global_dof - 1;    % 响应节点
    
    % 组装每个单元的矩阵
    for elem = 1:geometry.N

        % 获取单元矩阵
        D_local = [D_LL, D_LI;...
            D_IL, D_II_La];
        
        % 计算全局自由度索引
        start_dof = 2*(elem-1)+1;
        end_dof = start_dof + 3;
        
        % 组装到全局矩阵
        D_global(start_dof:end_dof, start_dof:end_dof) = ...
        D_global(start_dof:end_dof, start_dof:end_dof) + D_local;
    end


    % 已知条件：左端点位移 = 1
    known_disp_dof = excitation_dof; % 位移已知的自由度
    unknown_dofs = setdiff(1:global_dof, known_disp_dof); % 其他自由度
    
    % 分区处理系统方程
    D_kk = D_global(known_disp_dof, known_disp_dof);
    D_ku = D_global(known_disp_dof, unknown_dofs);
    D_uk = D_global(unknown_dofs, known_disp_dof);
    D_uu = D_global(unknown_dofs, unknown_dofs);


    % 已知位移条件
    U_k = 1; % 左端点单位位移激励

    % 求解未知位移
    % K_uu * U_u = -K_uk * U_k
    U_u = D_uu \ (-D_uk * U_k);
    
    % 组装完整的位移向量
    U_full = zeros(global_dof, 1);
    U_full(known_disp_dof) = U_k;
    U_full(unknown_dofs) = U_u;

    % 提取输入和输出位移
    U_in = U_full(excitation_dof);
    U_out = U_full(response_dof);

    T_N = abs(U_out/U_in);


    % 
    % % 计算透射率 T_N = |B11 + B12 * A21 * A11^{-1}|
    % % T_N = abs(B11 + B12 * A21 / A11);
    % T_N = abs(alpha_N(3,1)/ alpha_N(1,1));

end

function kb = calculate_wave_number(omega, material, geometry)
    % 计算波数kb
    % 输入参数：
    %   omega : 角频率 [rad/s]
    %   material : 材料参数结构体
    %   geometry : 几何参数结构体
    %
    % 输出参数：
    %   kb : 波数 [rad/m]
    
    kb = (material.aluminum.rho * geometry.A * omega^2 / ...
         (material.aluminum.E * geometry.I))^(1/4);
end

% 梁谱元矩阵La、Lb %
function [D_beam,kb] = spectral_element_matrix(omega, material, geometry, L)
    % 输入参数：
    %   omega : 角频率 [rad/s]
    %   material : 材料参数结构体
    %   geometry : 几何参数结构体
    %
    % 输出参数：
    %   D_beam : 4×4梁刚度矩阵
    
    % 统一调用kb计算函数
    kb = calculate_wave_number(omega, material, geometry);
    kbL = kb * L;
    % 计算分母Δ
    Delta = 1 - cos(kbL)*cosh(kbL);

    % 计算各系数
    alpha = (cos(kbL)*sinh(kbL) + sin(kbL)*cosh(kbL)) * (kbL)^3 / Delta;
    alpha_bar = (sin(kbL) + sinh(kbL)) * (kbL)^3 / Delta;
    beta = (-cos(kbL)*sinh(kbL) + sin(kbL)*cosh(kbL)) * kbL / Delta;
    beta_bar = (-sin(kbL) + sinh(kbL)) * kbL / Delta;
    gamma = (-cos(kbL) + cosh(kbL)) * (kbL)^2 / Delta;
    gamma_bar = sin(kbL)*sinh(kbL) * (kbL)^2 / Delta;

    % 构建刚度矩阵
    EI = material.aluminum.E * geometry.I;  % 弯曲刚度
    % L = geometry.L;                         % 梁长度

    D11 = (EI/L^3) * [alpha, gamma_bar*L; gamma_bar*L, beta*L^2];
    D12 = (EI/L^3) * [-alpha_bar, gamma*L; -gamma*L, beta_bar*L^2];
    D21 = D12';  % D21是D12的转置
    D22 = (EI/L^3) * [alpha, -gamma_bar*L; -gamma_bar*L, beta*L^2];

    D_beam = [D11, D12; D21, D22];
end

function D_0 = calculate_DVA_stiffness(omega, material)
    % 输入参数：
    %   omega : 频率
    %   geometry : 包含材料参数的结构体
    %
    % 输出参数：
    %   D_0 : DVA动刚度

    % 读取参数
    k = material.DVA.k;
    m = material.DVA.m;

    % 计算：
    D_0 = k - k^2 / (k - omega^2 * m);
end


function plot_k(fre, k_record)
    figure;
    plot(fre, k_record, 'b-', 'LineWidth', 2);

    xlabel('Frequency (Hz)', 'FontSize', 12);
    ylabel('k', 'FontSize', 12);
    title(sprintf('k_record'), ...
        'FontSize', 14);
    grid on;
    xlim([0 max(fre)]);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":26.1}
%---
%[output:599fbcb9]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAARkAAACpCAYAAAAfrUnWAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQu8llP2x1cXXVAUUVONihBFEaJGRSK3MGMoxpSiGlO5VlK6UJqUQUw3XfhLMySmmKJIuQ2JMhKGLpMkNCQlKuc\/333sc57z9r7nfZ\/nfW7v29qfz\/mcc9732XuvtfZev2ettfdeu0xBQUGBaFEJqARUAgFJoIyCTDCS3bVrl5QvXz5p4z\/99JN888038uOPP0rNmjVLPLNq1Srzf6NGjUp8\/u2338onn3wixxxzjFSsWHGPdr\/77jt5+eWXpXnz5lKjRg3XTK1YsUJq164tBx988B51V69eLW+\/\/bZceOGFUqFCBVdtw8\/mzZulZcuWUqZMmbR14eP555+XX\/3qV2n52LBhg\/B8w4YN5d\/\/\/rdp+6ijjiq1j6+++kqod+yxx5oxsH8nG6vPPvtMKlWqJNWrVzdt0tfTTz8tbdu2lX333VcOPPDAjHhKy3SeP6AgE8AAL126VEaOHCl\/+tOf5LDDDpOxY8cKnwEqH3\/8sdSpU8coRqtWrYziVq1aVQClr7\/+2tTh89atW0uVKlWKgAolok3aYtJ\/+OGH8uWXXxZRv3PnTpk+fbppu0OHDiW4Ovroo6Vs2bLSrVs3AUx+8YtfyO23324Uf\/fu3QKA0e9VV10lhx9+uFGcAw44wCjVp59+Kv\/973\/liSeekFGjRsmOHTtMXT53Fj6j7vvvv294seWjjz6SGTNmSN++fYuUle\/gDUWnPeRiC\/RAC2Dx61\/\/uuhzaHLKgz6QR7Vq1eSPf\/yjjB492jzbv3\/\/Ukf0n\/\/8p8ycOdPwgizs35UrVy5RDwP\/gQceMJ\/RPv0jh5tuusnI8bHHHjP89u7d2zXwBjDlYt2kgkwAw8MEnTVrlkydOlXuu+8+80bms++\/\/15GjBghl1xyiZxwwgmmZyyD\/fff37yJr7vuOmGyo0Dbt2+XSZMmGTCiJILM\/PnzxVo9fM8b+rXXXpOzzjrLWDqAii0dO3Y0b91hw4aZPsaNGydHHHGE6Wfbtm2yYMECo+wAFIXfw4cPN20+\/PDDcsYZZ0j37t3Ndyj22Wefberdcccd0rVrVznyyCMNDyjh7Nmz5YcffjDPYrG98sorUqtWLaOQFEtX3bp1DcAio+eee858B9jRZ7169Uxb\/M\/f1AGIhw4dKvXr1zfPvvPOO4bGP\/\/5z+YZgMmCDLJ+4YUX5NBDD5UmTZqUGOFMQWbt2rUyZMgQIzPap1ier7zySiNPaOBlAG1aUktAQSag2YECr1y50ry9586dawCGNydKd9BBBxkFQzFxi8455xwDIlgpv\/3tb42Sony8yf\/2t78ZCgEBrIdf\/vKXcueddxoL55FHHjFvUvrgs\/POO88oNnUGDRokjz76qLGkLr30UuMaoDAA3Ouvvy49evQwfVx++eXyl7\/8xbhh1gIC6PiBJgsyf\/\/7383bv1y5cjJnzhz517\/+JUuWLDGghssH3fvss4+MHz\/euF3wh\/X2zDPPyPXXXy9jxowxYAXgTJs2zVgcPEfZuHGjoRnrBfCkH\/pYtGiRnHnmmYIl5gSLL774Qvr06SPHHXec3HrrrcbKcIIMMr7hhhsMoEPfG2+8IZ07d0470hdccIHpGz54OeA68jKAD9yjp556Sj744AMzloAmlhigbYEvbQd76QMKMiEMPG9AFGv58uWydetW84bFRUDReEuiJCj04MGDzaQGoF599VXzv3U9iB1govfs2VNOOukkA1S8rbEcTjnlFGNhED+gbRTkH\/\/4h3GFACEAA8sAkMEyQikAJGtx3XbbbQYUUHoUFMBo3769oQ8FA8SshQEo4kYRMwEkATVAFOXkOVw+lBsAWLhwoQEjABK627RpI\/fcc4\/ceOONcvzxxxvJwx+WFbQ3a9asyH0BZLDMcOHuvvtuQwvf0zdgDPgBlNY9siADH7g3v\/\/97w2wAMZYLwAyfGL9ITP+xuW04Ol0nd588025\/\/77TewFGn7zm98UxXDoF5patGgRwszJjy4UZAIYRxR88uTJxu3ASkHpL7vsMvPDmxogAXiwPlDGm2++Wf7zn\/\/sATIoLhYIysXb3hmTgWxA5N133zWmO889+OCDRnkBCz5r0KBBkalvQQbL4L333jMgBwAAOps2bTJxICwJLAr+vuiii0zbWFVYLJTTTz\/dKC2KiQVEnOb88883wIBrYQOkW7ZsMYAK3Vg20IkFhCzoCxfDBq+XLVtmLIQBAwYYumyMxIIMMuCZefPmycCBAw3NuElNmzY1cSMnyLz00ksGhAAB3DgbzE3lIjmBxfk3cSRcSH7jCuEWQi8yBuCwpOAHunr16iVXX311ALMof5pUkAlgLNevX28UBuXhrc6qEBMVIAF4sAqIL+CeYAmgDFgyU6ZMMZYHAIDS4goAWNQBmCzIMPFRKKwCgAF3hbcz1gUrQShyp06dTF2U7eKLLzbxDSwZXBcCmigiCk0Qk7q4Bk6QARCJS+BS0T8KDXihePSBC5YIMoAJbWOFAVL8xorBcmLlC\/DCqkJJoYX2sKaw7rBokAsuIa4RlhvAiBu33377maCvtXyIVwE60AfIAIQAD64Z7WLNOGNSgAxAf80118iaNWuMhZXsb0AVUING+ARQaA+Aufbaa+Xzzz83NABwxKiQGcCaGDQOYErldJMKMgENn7Uc7BsehQAkUAZWW7AO+GzixInG9QEMWDK1hRUgXCyUj+AmYDBhwoSi1SVWcZjcuFvECVBQgIQgLEvGfIcS4O4AZlhOKEy\/fv2MdWFBArcCkGHZGPBCybGsABlcDApKBoCgoAARVgrxIpZzAStcC\/gkiI2iYsHwPM\/ZlSPosoFSlB66+QyAwNV58sknTf\/OJXS7GkcfgIlzmRmLC3rgDbnCp9N9cg4r\/bH61a5dOwNw0Jjsb1boABmsRlaRkAVuHK4l\/xOrefHFF42soIcXyGmnnWbiXFpSS0BBJqDZ4QQZFB6gYLmV2AQxGJQeawYlJ7DJ2xzFxPrhTYuCsSyLFcHqyyGHHGJcMN6uxEJwJ7CYbEFBqIf7ZQOq9ju7uoSiEazEqkBBsJwsyJx44onGWiI+Aj3WzQPciNMQp8BdAExY9Tr33HPNChqBTwsyfI8r5SxYJFapUVJb7OqS3XdDH1hBgLAFH0CWIDUxGfh3FkAGoGXljjgUliPBWes+WYuFlaBUblGqz7E4CToT68KdRN62YHHituEu8T2WZ6r9UAFNrZxrVkEmoCEDZIgzEEQlVoICs+KBImC9YNqjOLg8gA8WC4FdYjYEfAEYJnOyJWyvIIMlg4XCihPKzf+s0gBsxElYNYFmPgdkbEFh+RyFIuCJK7Fu3TpDa5cuXcxKCzRhmXkFGfrCGsFFxDogKAsg4t4lWyIGZKy7RF3n6hL\/QyMyBKBoK9nemFQgw0sBkMfNBeCxthg7wAQ3FRp\/97vfGeuSF4MThAKaTjndrIJMAMOH747VQBCSmATuA8qduJsW94YgKW9j3uDEHlh6Jk6BeZ4KZGyA1Um6BQhiMclWPpyWFYB2yy23mFUqVk8AFOI5do8OrherSLytcetYfibASZCXgmWCRQXAEPdhzw3uXrKdyM6ga7LYhd39DCgsXrzYgC\/PsUqF3LDCiIEgI+du40SQwQoCTJAzBVcNdwe5slsZVwoXjTGxK2WJfzdu3Ni4S7hEgCr\/YwkRIyLg\/NBDDxkaAD9cW+I6zz77rJGlxmVCcpeSTXT7GcucrHxgolslYaLwVqDwFs+XZUFMfyY5k5DC25A3Kbt9iZfwRuTtj7\/PxjaCiOypwU0hbkEMguMBBCt5BsDhDUocBRBKNqHdgAxvafoBXAiiEg8hkAogUlBqLAhcPBQRa4f9ObhkxIXgDZoIyvKW5zOAhthFYikNZAAYXC1iVLhr7FM5+eSTDaDAz1tvvWVcQOQJ4GEB2j0piSBDMBqAsLugkRvxJ4DKuYSdCgwAIrsUD10AvbVQADxcRJbgAWZkhmyYs3ZzZQDvqrxp0jdLhjel3bbuBAznZMCk5W3C28a5g5W3WKrt3bkoaSYoP84VjlzkQ2lWCfghAV9AhrcObyQ2LbEyQfANq8S+XfGxARbnW825D4LnEt0DP5jTNlQCKoHoJeALyFg2rDWTCDI2TgDIYM3gMrHZyQbuEutFLxalQCWgEvBLArEGGQJ3iad9\/WJc21EJ7K0SYFuEPQwbhgxCARkv7hLgQtSew21aVAIqAf8kwP4elvbDAppAQQaxeA384lpxwA1hJG4u80\/c4bYEYLJilE88IUHlK9x5lE1vdqzCXM0NHGScS9j2KL1dRrRL2HYLvd14hhAtyIQpjGwGL5O6xKA4xcsuW5ujJJN6cX9G+Yr7CBXTF4Ve+Qoyfoo6CmH4SX+yttiHwrkYcqqwjyNfivKVOyMZhV4pyIQ4P1QZQxS2D13l43gpyDgmRhTC8GFeltpEPk5aGFa+gp45\/rUfhV6pJePf+KVtSZUxrYhi9UA+jpeCjFoysVKyTInJR2XMVwtNQUZBJlO9jtVzCjKxGo5SiVGQUZDJndnqoFRBJneGTUFGQSZ3ZquCTE6OlYKMgkxOTly1ZHJn2BRkFGRyZ7aqJZOTY6UgoyCTkxNXLZncGTYFGQWZ3JmtasmkHKtu3USmTInnUOYlyDjTciYehCwtx28Uwgh6WugbP2gJ+9u+1\/EqU6aQjoICf+nxo7Uo9CrwHb\/OqyoAFbLjk\/CZfDH2RsRkOX6jEIYfg1haG14nbdB0Zdu+8lUswWHD5H\/3aRX+v2aNSL162UrX3\/pR6FXgIOPMJ5P4twWcZDl+oxCGv8O5Z2uqjEFL2N\/2vYyXgsyeYxA4yNCldYvuuuuuontxnICTLMevgoy\/ChNka16UMUh6\/GrbC18KMiGDTOJtBU53idsKSkskbkGmb9++5oIvLpXP9cKk5eI3rlXNp8vAlK\/imXnnneVlxIjy5oNVq3bEyl1i7i1btsyktQ0zGVyglgwWivOqE+5a4sItrg\/lMq1M3CUGi0xyXHiW64WL0Lh8jMvLkt22mKv8KV\/FI3fffQf+L8VqNfPBkiXrpU6dXbEZ1kceecRkZqSEDjLOFaB0Ekm8BbK055Pdu2SvROEmxUwCv+TD5XbBfLBkkMemTZsML\/mUGU\/5yh1LhrvKyTMdCchwyfqQIUOKrpBNBh7O+5ST3cecqo69WTLVEvbekuPXi4+fDvTj8L3yVTwKGpMJOSaTjQJo4Dcb6YVbV0FGQaa0GVcUk3G6TLhEY8aMMVfOrlixQpJZGkFPYwWZoCXsX\/sKMgoyaUEmMXZC3GTu3Lkybdo04ZoS56pQWKsiCjL+gUDQLSnIKMikBZnEWAurQEShb7vtNrPU6iUWk+3EVpDJVoLh1VeQUZBRkAlP30rtSZUxJgORIRlexssZ+F20SKRNmww7C+mxKF7eJiajlkw4I+xl0oZDWXa9KF\/JLRkFmUK5FIGMXWZONd3c7I\/JbsoW1o4Ccf2gu7Q2VBmDlrC\/7XsZL7Vk9hyDQHf8ZjPkCjLZSC\/cul6UMVwKvfXmhS8FGQUZb7PNp1peJq1PXQfajPKl7lJpE0zdpUDVr2TjqowhCtuHrryMl1oyGVgy7JGpV69eUUoGGx9ZvHix9O\/f34ehy6wJdZcyk1McnvKijHGgOx0NXvhSkEkDMqn2w+g+mXTTMbPvvUzazFqO9inlK7m7pJnxCuWyR+AXS2bixIliE0zZ\/3v06KGWTJa6rMqYpQBDru5lvNSSycBdsu5R586dzdN+nFuyQEV7ziPmmkg8ZK0JqDsvyhgQKb4264UvBZkMQcbPkXKm2XQeV9BE4pX8FHOkbXlRxkgJzrBzL3wpyKQAmUxjLpk+Z7vh4OWIESNMZjsOWjqL89ClJhLPcNbH9DEvyhhTVkqQ5YUvzSdTCsik2\/Frq7rZ+WtBxrpJTndJE4mrJRN3oPECMhMmiPTqVciZHisolEOgO35tjprLL7\/cLImzLG3Tby5YsEATicddyzKkTxOJFwvKmUh83rwdsTogmZeJxG2emk6dOkmLFi1EE4lrIvEMcSsWj3lJkO5MJP7YYxulRYsdseAFIiJNJB6kFJxukdOS0UTi+eMuaSLxYg3q3r28zJhReCXKggW7pFWr+NxWgCUTWSLxIEHGWjNk2tNE4jtk48aNUqtWrby6rcBL7CLIOedX21746tpVZPp0jck4xyDQmEw2g63HCrKRXrh1vShjuBR6680LXwoyKVaXvA1BsLUUZIKVr5+te1FGP\/sPqi0vfF17rcjkyWrJqCUT1KxM066XSRsRqa66Vb6KxaWWjEtL5rvvvpMff\/yxRK0yZcrIAQccIGXLlnU1Ed0+rJaMW4lF97yCjIJMabOv1JjM\/PnzZdWqVUX1P\/nkE\/nqq69k\/PjxUq1a4X2\/QRUFmaAk63+7CjLFMu3TR2TcOHWXXLtL7BeYNWuWWRnp1auX7Lfffv7P1IQWFWQCF7FvHSjIqCXj2ZKh4po1a+Tee++VLl26SNOmTQV3KYyiIBOGlP3pQ0FGQcYzyMyePdvcHtm1a1fZf\/\/9TTsVK1aU4447zvwOsijIBCldf9tWkCmWZ9u2Ii+9pO5Sxu7Shx9+KF9++WWJGakg411BVRm9yy6Kml7GS0Fmz5HSzXghzl4vkzZE8jx3pXypJZPWXXKbJ8bzbHRRUd0lF8KK+FEFGQUZBZmIldB2r8oYk4HIkAwv41W\/vsjatRqT2SMmE4YlY3PLcK0KaR8omuM3w9ke88e8KGPMWTLkeeHLCTKzZ4tcfHG8OI3CQwjtcjebTNwmEie3zMiRI2Xs2LHy0UcfycyZM2XUqFFSuXJlMypRCCPo6eBl0gZNkx\/tK1\/FUnSCzIgRIgMH+iFh\/9qIQq+KQGbYsGEyZMgQqV69un8c\/dwSjD3zzDOydetWsQmsNMdv\/uSTUZApVhnnNjLOMU2d6rs6ZdVgXoKMdcX69esno0ePLgEya9euNXc5JXOlohBGVqOXQWVVxgyEFKNHvIyXE2S6dBGZNi1GDEXkIQRuyeAmtW7dWkhAPmDAANcg07dvX+nYsaPUrFkzXqPlgRrNhetBaBFW8TJelSsXW6hkxSM7XlxKXub4tRbKihUrSsiZuMy6devktddeM3GY0q5EoSJXqlx11VVxGSvPdHjJGeu5sxArKl\/Fwm7QoH7RP3Xq7JIlS9aHOBKld5W3OX4t28mSimcS+L377rvlxBNPzAtLRnPhxkbfMiLE7XixdN2oUclY2\/ffxyeReN7m+E0FMnxul7CTXYWrMZmM9CAWD3mJXcSC8DREuOULkGF1yVkKCuLFaRR6pccKQpwDbidtiKRl1ZXyVSi+ZCCzZo1IvXpZidfXygoyDnFGIQxfRzNJY6qMQUvY3\/bdjhdWS2LCSAWZgG+QzGbIFWSykV64dd0qY7jUee\/NLV+33y5yxx0l+1OQ+RlkyOX7yiuvSNu2baVChQomMRWfkebBpt8899xzvY+Wh5oKMh6EFlEVt8oYEZmuu3XLl03zgHtkzy8pyPwMMiw1T548WapWrSpNmjSRli1byuOPP26Shf\/0009mcLjLOsyiIBOmtLPry60yZtdbeLXd8mU34rVpU5y4avhwkcGDw6M5XU9R6FXRZjxA5uqrr5YJEybISSedJCtXrhQ2wj355JMKMulGLsPv3U7aDJuN\/DHlq3AILMhwZom7l7BmAJxFiyIfoiICIgeZiy66yFyhunv3bpkzZ465DqVKlSrGfVJLJvuJosqYvQzDbMHteFmQIQDsvH8pTsvYkYEMsZdJkybJMcccY05E161bV5YvXy7XXHONvP3222rJ+DSz3U5an7oNvBnlq6QlA6h88onIEUcUfh6nuEwkILNx48aCm2++WbZv3y4jRowwYIM1w1mhLVu2mNgM5eCDD5Y6derI8OHDA79zif6iEEbQ2qjKGLSE\/W3fzXg5g76AinPPzF5vyezcubOA+5TGjBljRohgL6tM06dPlxtvvFE2bNigloxPc9fNpPWpy1CaUb6K4zHOGIx1n+J0GjuKl3eJwG\/Pnj1N3pc2bdrI0qVLpaCgoOiaWo3JZK+vqozZyzDMFtyMlzMeY2l03lwQF2smcpD5wx\/+IM8995wcdthhJibDfUu6uuTftHYzaf3rNfiW9na+OnYUmTOnUM6JYJIMfIIfkdQ9RA4yBH5xj7p37y733XefsAHv3XffVXfJp1mxtyujT2IMrZlMx8sCCZvwiMc4i3PvTByWsiMDma+\/\/tqciAZUSL\/JbZHs9OVyN3KFULy6S5Yp2kg8ba2JxEPTl0A7ylQZAyUigMYz5as0a+XWW0VGjUpu5QRActomIwOZtJR5fICdxDfddJMMHDhQGjZsaIDsr3\/9q0yZMkU2b96sicQ9yjVu1TJVxrjRnY6eTPgqzYqx7dtnLrpI5Kmn0vUa7PeRg8w777wjrDQFdU7JeUPBggULMsqMZ283CFb04bSeyaQNhxJ\/e9lb+dq6VaRq1fRWCsf+5s1L\/5y\/o5K8tUhAZuvWrQXshWFPDPEYfhObSVZOOeUUOfXUU80OYC8FS8YmD3f+rYnEvUgzPnX2VpDJxIpJtGaSxW3CHMlIQGb37t0FAAvL1ezuBWguuOCCPfhmsx7A0KNHDxOzcVtgznm3UqYgo4nE3Uo6\/OcBGVI7HnrooUX3ZoVPhf89lsZXhw6V5KWXCvvMJMXmG29wjqkwNWfv3rtk9OjwE4zHIpE4QECCb4K83JH08MMPS8WKFc1qEwJ\/8MEHpVu3bq53\/DoBxU6FTO9d4nlNJO6\/AvnZ4t6WSPyJJ6pI\/\/4HGxG2abNdpk7dlJE4O3WqJW+8UQg0q1cnLENl1EJ2D8UikbgFmWbNmpkdwL169TJXmbALGEuHH\/52UwATSuLqVKY3SGoicTfSjuZZtwm3o6HSfa\/J+Nq1S6RKlUKgwPVZtcpdonDnlSmZWEDuqU5dIxaJxFevXm1WfU444QRDably5bLiESBhQ99nn31W1A6gxeoSS+WaSDwr8cam8t4Sk2EPTIMGhWLPJrbiDGmGvRM4kphMAeaJiNi7qlPN3PLly8uVV15pbnzEhQq6RCGMoHnaW5QxaDmG1b5zvHBzOJeULcBY2p1A8\/LLIq1ahcNVFHpVZseOHQWbNm2ShQsXmpseScN51FFHGY4XL14sX3zxhQwaNEgqVaok999\/v0n\/EMR92YkijkIYQQ+zgkzQEva3fTtep59eVz79tHwRwLBz148bCJxAg4VEeoigSxR6VWbLli0FBIQ4EHn77beXAJm33nrL7NJ97733pHfv3ubO6vr167uOy3gRXBTC8EKnmzoKMm6kFf2zjJczhpKNi5SKm8TdIEG7T1HoVdG9S+PGjTOJw1ldqlGjhpHJN998Y9wjLB0OTHJKG7cpjBKFMILmS0EmaAn71363biJTpxa3FwTA2Na5gfn\/\/q9kX35ZS3HwEIpAhtsJSLfpLNOmTZMLL7zQWC8PPfSQtGrVSho3buzfSJbSkoJMKGL2pZN8As+vvhL5+R1rZMN91h99tKuEReOL0JI0kmjVBAFsUehVxjdIcmsBO3297vZ1OzBRCMMtjW6fzydldPKeD3xxDrhSyWusDcAsWbLeZIokJhlGWbJEpHXrkj35CTZR6JUBGQCEDVXPPvus2e1LwLdatWrmihRWkojLsBN4yZIlZmNcGCUKYQTNVz4oYzIZ5TJfuCm4K85ilTpKvvr0ERk3bk9pP\/usSDZXoEWhV2Z1iZPR5513nsydO9f8JhDMzt4ZM2bItm3bZM2aNXLFFVfIm2++aa5JCSMuE4UwFGS8SSBKZfRGsYgza51tI9FiiANf8+eLdOiwJ5fQunQpubfdSSAKvTKWDMcH2rVrZywZVpBOPvlkEwTmVHa\/fv1MUu\/99tvPrDyR97dy5cruOPPwdBTC8ECmqypxmLSuCM7w4Vzgi8Teo0eLjB+fXGETk03xVNz4ql+\/+GbKRMtryhSRM85IP2BR6JUBGSwW9uQBMk2bNjX7ZVi+ZucvRwzYjkwshh26ffr0UZBJP5ZJn4jbpPXIxh7V4soXwHLXXSKTJiUHlvbtRSZOTC2FuPJFssrjj09O99ChIkOGpOYpMpB5+umnTYqH+fPnS\/Pmzc01ta+\/\/rpceumlUrt2bQMw9ipbBRnvqhnXSeudo8KaceKLfSa80e0J6UTe3ARR48RXqjFavVrk8MOLv2Xp2+5MTlYnMpDBDWL5movdECy3FXCDJOkdAJfBgwebOAxL2mS602MF3tQyFyatF86i5itZfCXRnfCy7yRqvryMRbo6kYEMHQMobMQjNsO92NxYwFklPgNcLrnkEuNCcTLbr6I5fv2SZLTthKmMGzeKdO6c2lJBElgrWDPEKbIpYfKVDZ1u6kYGMqTc3Llzp1k9Ygn7hRdeMMvXbL6jzJs3z3zGChMxGj9KpqkeNP2mH9IOto2glJG9K+ecUzqgWFDhN66Dx6SNSQUUFF\/BjkbprUcGMpC1fv16QfG5PRLrhf\/txjsCv59++qlZvvZrM16mSasUZKKckpn17YcyPvCAyJNPpgcUCypYK0FfMeIHX5lJMLynIgMZ8siQ6uH666+XRo0amVUkjhnY09i4UpzS5s7sQw45xBeJZJp+U0HGF3EH2ogbZRwwQIRUlKkCs4mE2tPOjz4q0rJloGzs0bgbvsKlzHtvkYEMO3zZOn3kkUca6gEAYjItWrQo4obg8KJFi8y+GT8Cv5mCjOb49T6hwqqJMtocvxUqVJYhQyrKsmVlMgYSSyfb+GvX3iWzZ4tUrx4W9an7cfIVxt6woDmORY7foJl0tp+pu0QdzfEb5siU3tfWrWXlnnuqyQcfVJANG8oX5VlxQ6EFk0GDNsuxx5Y8lOumnaCfzbfcxbHI8Rv0oDnbzzTwqzl+wxkV9pcQE3nwwcL7nL0CiKUWN4ef9u0L5IYbCm8hzbWSb7mLY5HjN+xJoDl+g5X455+LTJ6g72leAAAJo0lEQVQs8uKLhf2wA5YfP4oFkUaNRIYO3SHbtm0M9bSyHzyka0NjMukklNn3Gad6yKw5\/56KIkDlH\/XJW\/IyaXfvFlm+vPB601dfLQYLCxp+02wDrfwmC+uNN4r8HKpL2ZUXvvymO4j28pGvKPRKQcbn2bljhwiH7TghyxmTZcuKO8CK2LVrl6c4RjZkOoGDdriTmT2VFSpk02px3XxURrjLR74UZBxzHmGcffbH0qZNC2nQ4AgTJ\/jpJxHe7PY3m6\/4O5O8qKW5CX65EP6obOpWnMmr+ZuNZ+xs7dKFDG5B9566\/XxURgUZ\/+ZTrC2ZU08tXkL3j+XwW7LggBXDEi3Jqbkj75RTCkGC\/R8hZM8IjHEFmcBE63vDaskkWDJuQSbTayqSPYfSE4OoUkWkbt3COETDhiK1au2ZltHryKsyepVcNPXycbwUZBJApnPnzqI7fqNRMDe95qMyqrvkZgaU\/mys3SUFGf8GOsiWFGSClK6\/baslo5aMvzMqpNYUZEIStA\/dKMgoyPgwjcJvQkEmfJl77VFBRkHG69yJtJ6CTKTid9W5goyCjKsJE5eHFWTiMhLp6VCQUZBJP0ti+ISCTAwHJQVJCjIKMrkzWx2UKsjkzrDlJchYphgGrrslKXlDdrn9nBzr1ltvNX8n7oeJQhhBTxVVxqAl7G\/7+TheUehVoPtkuE6FK1QGDhxogIXUDlyJO2XKFNm8ebOMHDlSxo4da65imTlzpowaNaro4rgohOHvFN2ztXyctHCpfAU9c\/xrPwq9ChRkEkXjTFS1YMECee211wywkBzICUbUi0IY\/g1l8pbWrl0rXAlMpr96mZ6BCJooH9pXvnwQYkhNRKFXoYKMM6+vmxy\/p3CSMA\/Khg0b5JZbbjG3PuQLTwyL8pU7k9OOVZjHdUIDGRDU6RKlAxmuYEEh3yC1vRaVgErANwnwgiOtbZ2Q8oP4CjK4Q127dpXPPvvMXBJnYyxOQLGSSpdInOcAGn60qARUAv5JAHAJC2Cg2leQSSYGwIRy2WWXlfg6XSJx\/0SqLakEVAJRSiBQkHFaNpbJ448\/3qwuVa9e3aw2sYSduLQdpUC0b5WASsBfCQQKMv6Sqq2pBFQCuSiBWIKMtXAQaJhRcD8GMNF669Gjh\/Tv3980nYov9hN169ZNVqxYUSKW5Qc9frSRGLQvjRfn5su77rqryE2OG4+JPLGNYsCAATJ37lwjMqfFnQs8edn0GhZfsQOZXI\/VJFNIJm1pfHEPOftmLrzwQjPRO3XqVOKKYD+AwmsbFhidgfxUvDj3O9Gf3WyJaxwnHpPxlLhx1MrL+XlcefKy6TXMsYodyGSy6uRVYcKol2wlzb75k20+POigg4wVg7XD3eOp6odBe2IfAOa6devMx5Z27oRONUbs4gZMiLnxnAVM7liPC4+peHICJ6BoC8\/HnafEcctk02uYYxVLkGEHKUpnTWyrgFEomps+E01uZ0A71b4gFDDx6IVTod30H9SzTlCxIJNsjOjf7oXib0DmtNNOk7POOit2PCbjyZ6jg3br6jkt07jzZMc\/k02vYY6VgkxQmvnz0Qj7FuQYRTLFVJCpHOAIpG46EWScTzpfbmEqox+CyHTTa5h8xRJkSjvT5MdAhNUGZmu\/fv1k9OjR8vbbbyc9qxVnd8n5ZswXdykVT845YS1SrLDDDjssZ9wlN5te92p3KZcDv0zOESNGmAOQ9tS5VU52Lqc6dR6noGgyAE586+d64DcxRmZdQD5n06jz5cBLwLqzfB\/nYLal3zmGcRir2FkydgLk6iY95xJ2qvw5iZ87l3edS95hWVzp+knmWqTaSOlcFnVuP4gbj4k8JcbTnMvvcefJ66bXsPiKJcikm\/T6vUpAJZA7ElCQyZ2xUkpVAjkpAQWZnBw2JVolkDsSUJDJnbHyhdLvvvtOnn76aWnXrp3UrFnTVZtff\/21vP7663L66afL\/vvv76quPrz3SkBBZi8Z+127dplcyvvuu68MGzasKO\/yhAkT5LnnnishhbPPPlu6d+8uK1eulK1btxZ9x98PPPCAdOjQQZo2bVr0efny5eWYY46RqlWr7iXSVDbdSEBBxo20cvhZu7x+\/vnnm30f9jDmoEGDpFy5cjJjxgwBLFjCrVChgvmZM2eOrF+\/vohrNhOyktGmTRtTp2zZsua7ihUryiWXXFLCMrIrF\/agoTOnM0vGlFTnhZxidq765Nph2RyeLr6SriDjqzjj25gTZCZOnFhkySxdulSWLFkiL7zwgrFGateuLR07dpQGDRrIwoULTfrTa6+9VgoKCmT48OFy3XXXybvvvivvv\/++sNw+btw4c2zgjDPOKMF8sp2nicclMgEZ2yjA2Lp169gcHI3vSMePMgWZ+I1JIBSNHz9eJk2aZJKFffHFFyb94iGHHCKHH364sWCaN29u+iXVwRVXXGGUGReL\/SQACvWmT59unmnWrJkMHjzYHIRs0qSJsX5ow1ncgozd+YyFZYsz3YKCTCDTIpRGFWRCEXP0nXz44Yfy5ZdfmpsFcF06d+5s4jNYJQcccEAJl8hpMfz444\/mSASpKHh29uzZcvTRR5t8Kzt37jSAU6NGjT0YdAsy9sI\/21DiFnkFmejnkFcKFGS8Si7H6v3www\/mDFXjxo1l2bJlJrA7a9YsqVatWomVIuIugAzg8dRTTwkW0JlnninUx31q27atARtK+\/bt5fnnnzcn5llxKlOmTJFUkoGM85SzfTBZ6tVkO4wVZHJswjnIVZDJ3bFzRfmqVauMe3PllVcacGnZsqW88sor0qhRI9lnn31k27ZtJiazePFiAzJcm\/HWW29J3bp15dBDD5Xly5fLt99+KyyB4yJ98MEHwglyLB2AiSVxGwiGMK+WjDN\/izOvi4KMq+GO1cMKMrEajmCI2b17t8nzgtUAaHA18Pbt201AF\/eJsmbNGmOxvPjiiwZksHgef\/xx2bJlSxFRBI9Z7uZ7LCBbUq0uJd6zlS7wmwpg6EdBJpi5EUarCjJhSDlGfWCx4DYdeOCB0rt3b9m4caMJBE+bNk3YH0PAlxgNFkqYIIOIbFqMxPiMgkyMJpAHUhRkPAgtF6uwBI1bw+lilqlJX8CuXYLB9957r7Fw+B8XaujQoUl3A7tZcnbrLpHfhGC0szjjNWrJ5OKsK6RZQSZ3x84V5QRo2d17ww03mHgM8RPcF\/Lf9OzZ0+ziZYMdILNo0SITzMUNcpZsQMYVsUkeVpDJVoLR1VeQiU72ed1z4o5fZxDXDeO649eNtOL5rIJMPMdFqVIJ5I0E\/h+fQx1ex7surgAAAABJRU5ErkJggg==","height":0,"width":0}}
%---
%[output:6f16dfd2]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAARkAAACpCAYAAAAfrUnWAAAAAXNSR0IArs4c6QAAHJFJREFUeF7tXXtsVkUWPwgr7SoaIYSH4DMlEf9A0ygFFV2IZkULq1Go3SwVuhUfYBNLgdb3srSlpa4tqAEttSRSqmsiwSUxohI1FddFJVpUGqF1KxSRuooGKiibM3U+ppf7mu+7d+7M\/c5NCF\/buTNzfuec33fmzGvAiRMnTgA9hAAhQAiEhMAAIpmQkKVqCQFCgCFAJEOGQAgQAqEiQCQTKrzqK9++fTvk5+fDhg0bICcnR30HFLTY3t4Oc+fOhYULF8Ls2bMVtEhNpIIAkUwq6Gn4LpGMhkpJ8y4RycTMAIhkYqbQGIhDJBMDJYoiEMnETKExEIdIJgZK9CKZ999\/H0pKSuDKK6+Exx57DM466yxbqY8cOQLNzc3wwgsvwN69e2HYsGEwb948mDNnDpxxxhmJd3DVQ1tbGzz99NPwxhtvAP6cnZ0N9957L1x11VVw2mmnsbI\/\/fQTbNy4keWHeH15eXmsvuHDh7MyPT09UFhYCLfddhv7vHr1ahg7diysWrUKLrnkEjh48CA888wz8PLLL8PRo0fhhhtugFtuuQUeeeQRyskYYrtEMoYoym83rZHMRx99BIsWLYIJEya4EsyPP\/4Ijz76KGzevBlmzJgB06dPh127dkFTUxNMmTIFHn\/8cTjzzDNZN959912477774KKLLoKCggLIyMiAl156Cd555x2oqqqCW2+9Fb755hvWLhIcJmknTZoEn3zyCatvzJgx8MQTT8CFF16YIJkffviB\/YwktG\/fPsjNzYVjx47B\/fffD19++SVrZ\/z48bBlyxZ47bXXoLe3F5YtW0aJX7+GEWE5IpkIwQ+jaZFkzj77bFiwYAFcdtll\/UjCrl0kl+LiYqisrIRZs2bBgAED+hEKRkIYgXz33XeMYJBwamtrYciQIawcktTDDz8MAwcOZFHGK6+8AjU1NSwyufbaaxNNfvHFFyxyufHGG2HJkiWA5II\/4\/8NDQ1wwQUXJMquX7+etfHUU0\/B1VdfzX5\/\/PhxFuXgP+wrzS6FYUXB1kkkEyyekdfGSeahhx5iwx789kdn5FGIXQcxKigvL4fPPvuMOfqoUaMSxXAItXTpUsD\/MfrA6WMkG4xYbrrpJlt5kXAeeOAB9jd8R2wbSQIjkA8++IC1NXjwYEYy5513HqszMzMzMdTCdvHB34vDtd27d7PoCKMcIpnITc6zA0QynhCZVYCTDO81DkHWrl0LF198saMgGJ3cddddMGjQICgqKoLTTz+9X9mWlhb46quv4LnnnmNDF691ODzPgut0MFqxPlgfRiKNjY0s74Mkc8UVV0BZWVkignKrg\/8Nh1ZEMvrbJ5GM\/jqS6iEnGRzyoANiXuTyyy9n0QOPEqwVcqfduXOnY1ujR49mpHDo0KGUSQaTwTiMEknGSkhEMlJq17owkYzW6pHvnJiTmThxIrz44ossV8ITslgjllmzZg2rHIc\/+HndunXw888\/nzI0sfYAE8lOwyWsY+vWrSzBvHLlSqnhkpVkrMM0cchFK37l7SLKN4hkokQ\/hLats0uHDx9m09cdHR1s2ISJVbu1NHZJVuwezhJh7gOnnJcvXw6\/\/PILS\/ziMEfM9XR3d7MpZRyeYdSEQyK3xO\/UqVMB80aYv8Hhkt3QCvtkrQOny\/H3ONtFid8QDCiEKolkQgA1yirtCASjjzvvvBNuvvlmNvODP69YsYIlXocOHcq6y6ewX3311cSU8\/79+xlZ7Nmzp98Mj3UKG9\/Hqemuri5GZDg8c5vCPuecc9hwady4cYkpbDuS4X3C6AhzRTSFHaVlJd82kUzy2Gn5ph3J4IwOTgPjvyeffJIRi5VkUBi+eA4Xvn3++edsevqaa66Bu+++Gy699NJEUta6GA\/fnTZtGouYxASzzGI8pyQx1oGRCw7Fvv\/+exbx4MxSfX09W1NDiV8tzbBfp4hk9NdR4D1EIrIjmcAbogoJATpPJj1tgEgmPfUeldQUyUSFfITtEslECH4aNk0kk4ZKJ5EJAZUIEMmoRJvaIgTSEAEimTRUOolMCKhEgEhGJdrUFiGQhgiETjKYZMSDkPgOW75cHI8WwAfPOREXhaWhDkhkQiDWCIRKMrhaFHfW4gFEnGRw4xsu2sKjBbKysmINLglHCBACId67hBFMZ2cnw7i1tTVBMri5raKigh1GxJe0kyIIAUIgvgiEGskgbBjNiCTDoxsOKW1yi69xkWSEACKgnGRE2PmZIXiwkXgRGW60w3\/0EAKEQPAI4BnL+E\/VEynJ8CTw5MmTExvdkFxKS0vZAdT0EAKEQPAI4DlDeISGKqJRTjI4XMIHd89ifmbx4sVQXV2dSALzXcQIwrnnnhs8whHXiORZV1fHlBxH+RBekjFiI3NpnutG5TXGyknGOoVtzcnE\/XIyPDwKz17BKz7Ek\/n1NUv5npGM8pipeiMK\/wqdZGTBiwIE2T6mUh4vKMPDoPBGALyvKI4PyaivVqPwLyIZxfZADqgY8JCaM1WPRDK\/HXLtdeVGSHajpFpTjVMGHJJRBi01ZRcuBFi9+mRb7723vd+Mbpi9oEgmTHRt6iYHVAx4SM3prsc\/\/AFg2zZn4Ylk8vPZJe3i2pmQbEV5tbobZxCAkIxBoChXhxepYG14A\/DcuW3w\/PO5Sv2LIhk5XaZcmhwwZQi1qCBKPR47BnDDDe6RCr9S\/D\/\/ARg27CRklJOhnIwWDpRqJ6J0wFT77vd9lTL+8APAzJnepILE8uabAAMGOEtBJEMk49fGtS6n0gGjAiJMGb\/9FuD22\/2RyltvySFAJEMkI2cxmpYO0wF1ETlIGf\/3P4BbbvEmlQsv7ItUUnmIZIhkUrEfbd4N0gG1EcrSkVRk7O0F+OMfvUll2jSA554LFgEiGSKZYC0qotpSccCIuizdrIyMJ04ATJ3qTSpYpqFBuitSLxDJEMlIGYyuhWUcUFcZvPrlJaPXlDImaXH4s25d39SyqodIhkhGla2F2o6XA4bauKLKRRm7uzNg2bI+wnB6kEjwX2OjWlKx9odIhkhGkYuE20w6kExV1a9QVnaa9qRCJGOjoiiYNlyX6197OjhgHGXcvx8gP987r9LWBvD736u0KLm2ovAvWvErp6OUS8fRAa2gxEVGr7zKmDHHYc2aX2H69NNTtgtVFRDJ0HBJla2F2o6JJNPRAfCPfwDU17vnVXJzAR54AGDkSDPPBSKSIZIJ1flVVW4KyXz\/PcCf\/uQ8BOLJWrtVtabISDkZysmo8nul7ejsgF5DICSW7dsBRoxwh0xnGd16TpEMRTJKySCsxnRywE8+Abj\/fvdoZcYMgLo6OTR0klGm50QyRDIy9qJt2SgdEHMrc+e6kwoCh0OgVBbBRSljKoonkiGSScV+tHlXtQPu2gVw333uxILJWjyCMqhHtYxB9ZtIhkgmKFuKtJ6wHRCjlcpKgLVr7cXkEUqq0YobiGHLGJYCiWSIZMKyLaX1huWAbklbJJa\/\/hXgz39ObRjkF6iwZPTbfrLliGSIZJK1Ha3eC8oBP\/gAYNYsAIxc7B4kljCjFYpkgjErWvEbDI6+awnKAX03GEHBZGVEMqmuBnjmGWdSCeOMlWQgSlbGZNoK8h2KZCiSCdKeIqtLxgH9zAbhKlxcaavTIyOjTv0mkiGS0ckek+6LlwP6IZa9e5NuXsmLXjIq6UQSjcSSZFCo5uZmqKqqgszMTAZLS0sLlJWVsc\/W+5WiACEJXSX9iqnGKSOwnYx79gAUFtpPM+ty1kqqMsq8H1XZKPwr1JwMJ5Pc3NwEybS3t0NFRQXU1tbC7t27TyGgKEBQqfB0IpkvvxwFRUUZtolbJJabbgJYtEjNbFDQOjZVj1H4V2gkg8J0dnYy3ba2tiZIBomH\/3zkyBEoKSmB8vJyyMrKYmWjACFoA3Srz1Tj9IvRP\/\/Zd52H3YPEUl4OcP31ZhKLKJOpeozCv0IjGa4QkVRwuIQ\/d3R0wJIlS6CnpwcKCwvZZ34lbRQg+HWgIMqZapxusj\/\/fN9Sfidi2bABYNKkINDTpw5T9RiFf2lLMsXFxTBz5kwYOXKkPpYVQE\/QOLu7u2HEiBGJHFUA1Sqvorn5dzBv3kDbdvGslXXrvoacnKFGy+gVkZqmR+zvjh07oLS0NF53YdtFMn6GS6jggoICmDNnjnIHCrPB3t5eOHjwIAwfPhwGDx4cZlOB171tWybMm2dP+nhK3Asv7IexY4+DyTL6Bc1EGdevXw9NTU1MROuEi1+5kymnPJLxm\/itqamB7Ozs2EUymIc6cOAAkysjIyMZnSl959NP8XbDQdDVNeiUdjHHsmXLUXa1h\/iYJmMygJooI0YymzZtgrq6uniTDCqUzzqNHj0aGhsbE0lf\/FsUY8ZkjCzZd0wYy7utY0Fi2boV4OKLnREwQcZk9cffM1XGKPwr9EhGVplRgCDbx1TK62qcXsSCs0bZ2f4k11VGf733V8pUGaPwLyIZfzYVWCmdjNOLWP72N4C\/\/EVedJ1klO+9vzdMlZFIhoZL\/iw8xVIVFQAPPnhqJTgUmj4d4KmnUmvAVAeUkdpUGYlkiGRk7FyqLEYteB6L3bEJ110X7PWppjqgDKCmykgkQyQjY+eeZb2GQ2GdxWKqA3oCKhQwVUYiGSIZGTt3LLtqVd8J\/dYHh0P\/+hfA+PGBNONYiakOKIOKqTISyRDJyNh5v7JOUQsSy9SpAA8\/rG6\/kKkOKAO+qTISyRDJyNg5K4un9D\/9tH3UEtZwyKuTpjqgl1zi302VkUiGSMaXnTslcTFqwVkjPFA7ysdUB5TBzFQZiWSIZFztHI9JwKtA7HItUUUtdh021QGJZGQQ8F+WFuP5xyqQksk4oN1VIBi14HoWXNei25OMjLrJ4NUfU2XUJpLBc14+\/vhjmIoZQ5sHN4fV19dDUVERDB061EsfUn+PAgSpDqZY2K9x9vQADBtmH7XQ+bcpKiGA1\/3qMYCmAq0iCv+yjWT4YVJ5eXkwe\/bsfkLyTk6YMAEaGhqIZCRNwMs4MZ\/S0NC\/Uoxa8vMBli+XbCyi4l4yRtStQJs1VUZtSAa1wTtTWVnJiIYTz86dO2H+\/PnsNLswnihACEMOpzqdjNNpSKR71GInp6kOKGMHpsoYhX+55mR4hzj4YUUvonKjAEHGuFItazVOPIvFutQfIxcTyYVjY6oDyujWVBmj8C\/PxK81opFRRDJlowAhmX4m+w4aZ2vrtzBt2ph+VZg2JHKT31QHlNGpqTJG4V+eJCMOnVQc2RcFCDLGlUpZuyloJJdduwB+u5Iqleq1eddUB5QB0FQZo\/CvBMmIORc\/YIc1dIoCBD\/yplIGLzVbt65\/DaYPiSiSOQr79++HUaNGGXGMKtdXFP7lK5JJxcFk340CBNk++i1vl8zFA7fffvu\/xhmnX5mxnKnf8ukgYxT+RSQjY1k+y7rNFJED+gRR82Km6pFIxvBtBX6moU01ThmfJxll0FJblkjGUJKxIxe8UdGah6GhhFqHCrM1U4mUSMYwkrFL6DqRCzd4U41TxmFJRhm01JYlkjGEZNauBZg\/P7nZInJAtU4VVmum6pFIRnOSOXr01PUsslPRphqnjLOSjDJoqS1LJKMxyViX\/8uSCw2X1DpT2K2ZSqREMhqSTFDkQiQTtturrZ9Ixj\/etE7GAau\/\/73v8G3+YOSyZw\/AgAH+wbUraapxykhNMsqgpbZsWkQyeODV0qVLYfPmzQxd6\/aEKECwqtlKJI8+CvDYY8EYAzlgMDhGXYupeozCv5RHMrhHqqSkBMrLyyErK+sUW4kCBN6JoIdGFMmYta9HhriIZPyjpZxk2tvboaKiAmpra21P1YuCZLZvB5g0qf\/QKKzzXEw1Tv8mRXuXZLBSXTYK\/1JOMi0tLVBWVpbAlp+8x3+hGgRr9PLQQwDLloWneiKZ8LBVWbOpelTtX6gT5SQjGgI\/XgKP8szJyWF\/4iAUFxfDzJkzYeTIkaHYzquv\/g5uv31gom7cHd3efjyUtsRK0Ti7u7thxIgRkBmnQ2QEIUnG0M0oqQbQ7nbs2AGlpaWg4mwo3slISYYngSdPnpw4sFw88rOgoADmzJmTFKBuL02ZMha6ugYlimzb1gXnnXcs8HbsKuzt7YWDBw\/C8OHDYfDgwUraVN0IyagacX\/trV+\/HpqamljhWJMMDpfwwcPJMT+zePFiqK6uTiSBOcnU1NRAdnZ2oJHMsWMAZ52VkdAITkt\/9tlRfxoKqBQS64EDB5hcGRkn+xJQ9VpUQzJqoYZTOoGRzKZNm6Curi7eJGOdwlaVk7HmXt5+G+Caa9Qbg6ljeRmkSEYZtNSWTbucjB28YYAgrntJdjtAUKZADhgUktHWY6oew\/AvL01EmpMJm2R++QVg0MnUC0RNMCivqcbpZUji30lGGbTUliWSCXDv0uOP91+li\/kYkXDUqvZka+SAUSEfbLum6pFIJiCSEfMvOkQv9C0frIPrUBuRjH8txG64pDPB0HDJv2HqXpJIxr+GYkUyYoJ32jSArVv9A6GqpKnGKYMPySiDltqyNFxKcriEd0ljBMOfTz8FuPRStcrz2xo5oF+k9C5nqh6JZJIkGTGCwY2NmIfR9THVOGXwJBll0FJblkhGkmSsEcyJE2oVlkxr5IDJoKbfO6bqkUhGgmREgsHI5a239I5guJuYapwybk4yyqCltiyRjATJmDREEs2IHFCtU4XVmql6JJLxSTKmEgwavKnGKeOsJKMMWmrLEsn4IBlxHcyvv6Z+sLdaFRPJqMY7rPZMJVIiGQ+SEe+cxnum8UpY0x5TjVMGZ5JRBi21ZYlkXEjm8GE8C6ZPIbptFZAxE3JAGbT0LWuqHolkXEiG52FMJhjKyehLGrI9I5Lxj5gR2wrEPIwJa2Hc4DfVOP2bFOWdZLBSXZYiGZtI5vPPAS65pE8V\/\/43wBVXqFZLsO0RyQSLZ1S1mapHIhkbkonLMIk7g6nGKePMJKMMWmrLEslYSObBB3PgzTf7lGD6MIlIRq0zhd2aqURKJGMhmUmT+u5iuu66vm0DcXhMNU4Z7ElGGbTUliWSEUjmm282wE8\/9ZFMXKIYlIUcUK1ThdWaqXokkhFIZu\/ePbGLYohkwnJ59fUSyfjHXMsp7Ly8e6Cz86PYRTFEMv4NU\/eSRDL+NaQlyeD2gaNHc4xe2eukAlON079J0ZBQBivVZWm49NtwiSd8v\/gCYNw41WoItz0imXDxVVW7qXokkrGQTJwSvtz4TTVOGeclGWXQUluWSEYgGdP3KNFwaT+MGjUKMjIy1HqRotZMJdK0IZmWlhYoKytj5rBhwwbIyembqsYHQcDhUpzWxoh2b6pxyvguySiDltqyaUEy7e3tUFFRAbW1tbB7925obm6GqqoqyMzM7Ecy5eWHYfnyIWo1oKA1ckAFICtowlQ9pgXJYBTT2trKiOXIkSNQUlIC5eXlkJWV1Y9k3ntve78IR4HdKGmio6MDmpqaoKCgAC7Q+e6WFNAgGVMAL+RX04Zk0AiXLFkCPT09UFhYyD7zIRMfLq1c+TJcddW5IUOuvvqvv\/4aSktLobi4GCZOnKi+AwpaJBkVgJxkE1w31jRFktX5ek35OhmMZNxIpqurC+65531oayv1JQAVIgQIATkE8MutpqYGxowZI\/dikqUjIRm34RLKgUSD\/+ghBAiB4BFAclFFMNh75STjlfgNHlKqkRAgBKJEQDnJoLB8Cnv06NHQ2NiYSPpGCQS1TQgQAuEgEAnJhCMK1UoIEAI6IqAVybgt0tMRPKc+4ZBw7ty5sG\/fPlZk\/vz5bAZNjOLws5jh5zNtO3fuhNzc3H5rh3STHWcAreubnHTHp0xRhsrKSpg9ezYTR2d5rfLhUoulS5fC5s2bWd8nTJgADQ0NMHToULZ4ND8\/n\/3eFPlU25M2JBOnXI2dE6Ji3WRcsWIFWzczY8YMZtB33HGHluuEOJmIROgkl7gOCuXnizDROXWV104+JETrei5OlPz3psinmmCwPW1IxmuRXhTgJNumOE0v1uEk47Bhw\/qtF3J6P9n+BPUekmdnZyerjs8Q4kptJ7kOHTrEyAS\/9bEcJ89x48ZpKa+TfCKJIkHyB8ubJF9QdiBbj1Yk47Z+RlawqMpbQ2sxue20RgidTvymFJ2Wb7eISh67dq39c5IL3+XDKvyMJDN58mS4\/vrrtZbXTj6+104cFokRq0nyqbYlIpmQERe\/7V5\/\/XXbhYhEMie3muhAqm4kL65SN5VEQzb5U6rXimS8FumpBieI9jDUXrx4MVRXV8OHH35ou2\/LlOESx8Pum95Od6YNl5zkE+2AR6oYkZ1\/\/vk0XPLhJNqQTFwSv2iEy5cvZxsgcdOn6JC4itlpB7quiVA\/w6U4JX5RXjsSxd\/jzJj4pYFfDpT49WYZbUiGKxfHvqYv0hOnsK2yOC1EFKd0xSlvbxWqL2E3nHCSS5zidZqy101eq3zWPJs4VW2ifKotRiuSUS08tUcIEALhI0AkEz7G1AIhkNYIEMmktfpJeEIgfASIZMLHmFogBNIaASIZxerHWaQ1a9bYtiruiVHcLeOaExcAYvJ10aJFtjv6edIWE\/B8\/5ibsE5bCIwDSKMOE8koVgaSDG6cFA9PV9wF45uzLvMPkmQQHKe9Z8YDF5EARDKKgSeSSQ1wcTEc39EdNMnwNnTdpJoagurfJpJRjLkfkuEh+5QpU2DZsmX91g2JRypg1+0OhBaHZLgGBXd3iytyrYe3Ow0pxHqs6334domioiJYsGBBAkVrf8T1P1iIDwnxM\/YjLy8vcfwDjyKchj74d3ExnHjDhd\/hkhU\/3nG79UziJlDFZhKr5ohkFKvTL8mgA+LDzy3Bz+ggq1atSuQe+KK\/hQsXJhwV60cC4O9xouBHMyCh+CEZaz3WqzT4z+KRD9Z37PqHMmzcuJH179lnnz1l6OiFj91CwFQiGU6CeFuGmLOxIzPFphKb5ohkFKvSy4mwO9zwxW95u99x4uFOi3uF8LCslStXJs6i4VEKluV3XXmRDCcHsR58X+w7Hq5ljR6s73ntJreSg5OMXEVOEZe46tZJnXariq3YiJsz7a7rUWwqsWmOSEaxKp1ml8Rw3ek+KrshgfiNK25IFM89EcnBTyQjRhtO56fg7Z\/8LBVeRiQZHBbh0Q5uszpWUnGLSJzI12uI5Ta7ZI0MRVOwy\/0oNpXYNEcko1iVMpGM9dI7fsyjtcucoHCXN49qRHKwHirlFck45S3EnIpfksHdyjxBawe12Lf6+nrXmTenSCeZ4ZLXTYqyU9+Kzcio5ohkFKsrFZJxS4jyb3RrdGEd5qQSyYhQiefkJBvJYH08EkNCxb5bE8Fim0GRjFMehiKZcJyBSCYcXB1rTZZknPIkorPL5GTs8j08+ekUGYhRB+Zk3IZLWJcfWXnEMGTIEGhra2Pn7vBZIyuIbjkZv7NLVtJ1OiSLcjLBOQaRTHBY+qrJj+M5Gbg1h2D3zW7Np\/ChD58Fwk5irgQfviCQ54nE5KjTTBFPBntFMkgydrNLdgRm7aPb6Xipzi655WFEBdLski9z9lWISMYXTMEVSoVksBfWfIl4tgnvpVgGyQWjhMOHDydIxbp2BevA85XxEadxrUlqcQ2MH5LB+qxt2Z0VZEdGToinsk4G1\/RgPgqjMLtHlM9rZiw4i4h\/TUQy8dexr2FLlDB4zSqpzpXQit9grYFIJlg8tazNT\/QUZcexf9Yoyq0\/TleUBCUD7V0KCsm+eohkgsVTy9p0JRmZXIwVWHEXdpCg0y7sINEkkgkeTaqRECAETkHg\/2\/fALn3sFYJAAAAAElFTkSuQmCC","height":0,"width":0}}
%---
