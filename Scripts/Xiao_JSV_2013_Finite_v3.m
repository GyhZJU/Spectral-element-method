% 目的：找到FRF求解方法
% 文献：Flexural wave propagation in beams with periodically attached vibration absorbers: Band-gap behavior and band formation mechanisms
% 
% 状态：失败，不是动刚度格式问题
% 备注：尝试更换动刚度矩阵格式。
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
geometry.L = 0.1;        % 梁总长度 [m]
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
omega_values = linspace(10, 14000, 10000);
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
figure; %[output:65a2327f]
plot(fre_values, T_N_results, 'b-', 'LineWidth', 2); %[output:65a2327f]
xlabel('频率 [Hz]'); %[output:65a2327f]
ylabel('透射率 T_N [dB]'); %[output:65a2327f]
title('透射率随频率变化曲线'); %[output:65a2327f]
grid on; %[output:65a2327f]
xlim([0 2000]); %[output:65a2327f]


plot_k(fre_values, k_record) %[output:4faf0d88]
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

% 梁谱元矩阵La、Lb %
function [D_beam,kb] = spectral_element_matrix(omega, material, geometry, L)
    % 输入参数：
    %   omega : 角频率 [rad/s]
    %   material : 材料参数结构体
    %   geometry : 几何参数结构体
    %
    % 输出参数：
    %   D_beam : 4×4梁刚度矩阵
    
    % kb计算函数
    kb = (material.aluminum.rho * geometry.A * omega^2 / ...
         (material.aluminum.E * geometry.I))^(1/4);

    [s11, s12, s13, s14, s22, s23, s24] = calculate_matrix_elements(kb, L);

    % 计算分母Δ
    delta = cos(kb*L) .* cosh(kb*L) - 1;
    % 弯曲刚度
    EI = material.aluminum.E * geometry.I;  

    % 构建对称矩阵
    D_beam = EI ./ delta .* [ ...
        s11,  s12,  s13,  s14; ...
        s12,  s22,  s23,  s24; ...
        s13,  s23,  s11, -s12; ...  % s33 = s11, s34 = -s12
        s14,  s24, -s12,  s22  ...  % s44 = s22
    ];

end

function [s11, s12, s13, s14, s22, s23, s24] = calculate_matrix_elements(kb, L)
    s11 = -kb.^3 .* (cos(kb*L).*sinh(kb*L) + sin(kb*L).*cosh(kb*L));
    s12 = -kb.^2 .* sin(kb*L) .* sinh(kb*L);
    s13 = kb.^3 .* (sin(kb*L) + sinh(kb*L));
    s14 = kb.^2 .* (cos(kb*L) - cosh(kb*L));
    s22 = kb .* (cos(kb*L).*sinh(kb*L) - sin(kb*L).*cosh(kb*L));
    s23 = kb.^2 .* (cosh(kb*L) - cos(kb*L));
    s24 = kb .* (sin(kb*L) - sinh(kb*L));
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
%   data: {"layout":"onright","rightPanelPercent":30.7}
%---
%[output:65a2327f]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAUEAAADCCAYAAADXXXDzAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQe0VcXVgPejV0WF0AULaFRERRFBlzUKiRCjKIJ\/QETArkiRYgNEmqhgRQMEYsOCFHEparAgf0TAILEbQHwoBCyhlwfvzzf+83Le5d73zq3nnDd71rrrwb1z5szsmf2dvWfPzMkrLCwsFE0qAZWASsBRCeQpBMtuz+\/Zs0fy8vKkQoUKKTWyoKBA9u7dK5UrV07p+t27d0ulSpXiXsuzd8uWLVKlSpWEeUq76bZt26Rq1apSrly50rIW\/b5r1y4pX768L5nQ\/lRl57dC+\/btM33Ex6YdO3YYudeoUcNvMZovDQkoBNMQXrYvRRnmzp0r55xzjtSpU0fWrVsnq1evNrfdunWrfPXVV7Jhwwb56KOP5Pzzz5frr7++mNLOnz9f3nnnHRk5cuR+IJs5c6asWbNGbrvttrjNAAAPPfSQyTNmzBgDm9j04YcfyqJFi8zXRx99tPzmN78puj91Hzp0qHTp0kXatGmz37X8PnjwYOnatav5fdOmTfL555\/vl49ya9euvd\/3X375pdxyyy3mHqeffrqvrgCAd9xxhxx22GFyzTXXFANPbAE2L\/enjtS1bdu2pj0lJfrk3nvvlQkTJphsvXr1MjKOJwN+X7ZsmTzzzDNy1113yQEHHGCuoW\/++c9\/SpMmTeSnn36Snj17SvXq1X21UTMlLwGFYPIyy9kVWDoTJ06Ujz\/+WCZNmiQo5uzZs6Vp06byr3\/9S1588UXze4MGDYpZDevXr5dZs2bJW2+9JQcffLAcc8wxBhSnnHJKUd1jIcj\/hwwZ4qttLVu2lClTpsgbb7whb775phx44IEGxEcddZQMHDhQnnrqKfnb3\/4ma9eulUaNGhlrDwW\/++67pV69ejJq1Ci59NJLZerUqVKxYkU5\/PDDDZieeOIJwXr0XnfnnXcaGACTFStW+Krf6NGj5aKLLpJx48bJ8uXLi10DbDdv3mxk6LUgL7jgAgNGm4A7MHv00Uelbt26+0GQvsFKjLWSk4EgDxruQfu6d+9eBGXbN0D+6aefliOOOELOPPNMX23XTMlLQCGYvMxyegUWH\/BAEYAbQCShyPn5+XLkkUcWuZNWkVFErL\/hw4cbQGHRASh+t1PA8+bNM1beTTfdZPKg1FhXhx56qHEX+ffYsWMFCHFf7n\/QQQcZQAGsWrVqyQsvvGDqgtuGhQrMPvjgA+nUqZMBNNC++OKLTR6AAYyBnheCxx9\/vBx77LHmN0D52WefyQMPPGDuzf1IuIxYRkCL73AVH374YZP\/t7\/9rXzzzTfSqlUrYzVRdx4KgBd3mzrwwcrkgcHD5McffzQypW2NGzc28uNj3U\/Kufnmm+W0004zYNy5c2cxCFLeiBEjTD2QMRb6nDlzTF25FusbWZP4\/owzzpCGDRua\/3sfRv\/4xz+MLIA1csHap1+fe+4583Bp3bq1\/P3vfzcPEuqr1mB2VE8hmB25ZrRUIEACiFhKgAY3ECVHwYYNG2ZcRqvIQLBPnz5SrVo18x2w7NGjh7EMY60pa9X9+9\/\/Nq5i8+bNjZUGRI877jhzX6wirKslS5YY5cQ9BoRYLIDg66+\/Nu44+XHP+e7aa6+VE088sUj5f\/3rX0v79u1l48aNRul\/\/\/vfGzcQSJ500kmCe3zDDTcUqx8wmzZtmrGUAMOMGTMMfIA3AO7YsaOxRk8++WQDQ+RAubivdi6SumNhdu7c2ZSFBfvDDz\/IoEGDzDTD4sWLDYgAOMlOA\/DgoH64sdZ1xx2+5JJLTDnU57HHHjPywurt27dvkXUJZONBkLpcdtllxqW2DzfqUL9+fdOuP\/zhD6YuWMtYmJdffrl5OPAAOuSQQ5Ka+8zoACzjhSkEI9TBWGQoHlbRhRdeKM8++6yBHUqOcmGB1axZ08wVxlqCzG1557PizQmimK+99pqxql555RVZunSpASgWFFYSEAVa1kKjDBQe6+3WW281FuXPP\/9srBaACUgAHFYndWduDAtv8uTJ+0md66k3c2+4+r\/61a8M6Pk0a9bM5MflJtgDxKgb86WAAhkAQ65p165dkVuJhYaLjysNVLi3F4L8n4fCyy+\/LA8++KCZOgB8AJ8E7LwQxDID8LQLa9XORQJB+oKHA\/JKNCfI\/XDD6QdkhyUK8L7\/\/ntzrbWaKY+pAuYVsdwpE2vQzhlGaMhGoqoKwRB3E8rQrVs3U0Pm2oANrhjWCO4r\/0cxcbWwjLBUzj33XGOZMZn+3XffFbUOxcYyBCJYHM8\/\/\/x+gRGAhhuLsqGsAPe6664z81JYJF988YX069fPuMwkFNkmoIeFisUCeIDS22+\/bcpiHg5o2yAMsMW6whojLy498MS1xD3HZeXvPffcUwRBrGECMQAVqxZYYSkCZNxK2nXCCSeYB4Gdp8NFB5x8BzRjIYhFikULTE899VQz9\/f4448bt5O5xCuuuKIYBN99911T9n333VcsGJMKBHF7eTDgOiNXHlrI5dNPPzXfvfTSS8bqp920i9+xOjVlXgIKwczLNKMlel0x4IXCvv7668YCw2IDPLhSgAkrkHkzrBusCBTHzgliCaLoWB\/k5zobHca6wfLAogQwlIM1g0sJbHHXsOhwZwEArh\/uJYpKwu2dPn26ASZgIdL55JNPmsl8LCmv1Wnn0\/gN66l\/\/\/4GsgRtUP54EGS+DPgSob7qqqtMXYEb9eY7rNQBAwYY0AG+8ePHmwcDUAH2uPrIEfAw74bLy9wdAAX2QNJGv4GOjWzbyDVlYqmSAGksjLzuMODCcrMPCh5I1orje+ppLXLcY+QE\/JE1wSyAjuULcOknphF4wDBtQb9oyrwEFIKZl2lGS\/RCEOVBgVFEJvSx+lAc5gbff\/9949JhicWLpGIJ4jITUACMWDTeJTJAFdigaMCCYACWJnN4zPcREAGkzPmxpg0rzFqCABIXEVBRBwtBr9sLOLEEUXIirli2KDlzeICK+UisMCxYYAwcrSWItcZ9gTGQwl3GTWQJ0HnnnWfkTb2BGgEQ6kkgBghiMRJRZ84UkNiIMOUAIeSA5YkrbZOVOXXDMqNtWJrI0AKM+\/OhPCDIlEHv3r2NvLinjTTz0MCaJzCFvHChKYMHFUuasNaRO33JfCBtoN5Al4gxrrJ11711zOggc7wwhWDIB4AXggQsFixYYObUmJPCiiEBM4CBkmLpkFBQXFsUCJf07LPPNsr37bffyvbt203U0UIQ6wULk9+8CQAQeQU0KKpNQBSgAiGuwXXEpcUSpA6JLEGsINxZ3HKCGVh9wJd6\/e53vzOWK\/NitJGybr\/9duMOE\/G10VcvqGIjr\/a32OVAdk0hc4dYfQAQIGG98lCIteyszKkLLityZtmPBRj3od5YvdQTkCU7J8g1ALtFixb7BT3oO6zZhQsXmgcODxZ1hbOnqArB7Mk2IyVbhcS9wpXEggAOKA9zbVgMuGtABMhgAWEtouRYiSgT84V2Qt5Wyuuikv+TTz4xFqU34ZKtXLnSWDjenR9YTyxtwcKiTri8lE8ZgAHrB2jg8qLEBDuuvvpqM8dmE\/fEYiLajOuOBYsFBLixQJlD5B5YeMwXxi6kpq5eK8tbb6LJWMreBKztnCRuNCAH1nbpijdv7ELuWGucvJSHvHn4EKhJFoKUwcMH+b766qvGQgZ2WLtY3bjIq1atkj\/96U+mriz50R0kGVGp\/QpRCGZHrhkpFbcUxWCujGUSuFosPWGi3LudC6Vmvglrjcl9lBtwWLfZG5W0FWOJB1Zeoh0j5AOULCFJtGPEusNYcriWABcYE2H+4x\/\/aOYKsULtmkbWA7K+EAsUFxOIYzkCIqCIlQgYCZgwV4Z1eeONN8bdrYIlVtpuDCw+AAVACdJg7RKN5VrqxRIb1gIiW2Brt67FQhDLjHpSXxvBxZJ+7733jAtPYMMGsErreOtSA3+sSGTE\/CpBLcDKfCdrQZn7pM48DHgAIUO+05R5CeQMgix\/8EYs7fo0liV4k998mRdF+ErEjWO+imAGbi6WFq4wc1AotzcBEubZUGxkiDXEfBKgAJBYFKwpxNVkvgpgAtgOHTokbHgyECRyy7wZygysuRZ4AzcS81m4u1h5zAUCc+qKVYnLyeJkrEHW+tEWQETbWXsIWGNTaRC0ARiW+vDQwCpl7R5zmTw0kBEBJpb9ADHqA3B5uMRCkHsTiMHKZa6URL4rr7zSrE8kyMQaRuoeb3uhrTvgt0uVCDhRD7vsBTkBQFx\/phqAMtYhlj7TFzwUCJJoyrwEcgZB7zKCkgaK33yZF0U4S0xnEz+WEC4mc3DZOAgAWJBK6s9YqXoDCt7f4rUz3uECyfQSUKbd3sMJ4l3PffhYGWXicIdk6ql5g5VAziBY2oZ9Kwa\/+YIVm95dJaASKCsSyBkEY3cK2C1JsYL0m6+sdIC2QyWgEghWAjmBYGx0DZeXyBrbkuyWKMTgNx95mefio0kloBKIvgRYUM4niJQVCHq3e8ULgMRbchCv8YnyAT+CAExIa1IJqASiLwF2M7G6IAgQZgWCpXVJvOhbSRC025dsHgtZuz2qtPuVhd8BPrsStM1loTcTt8Hlfk40RZbtHs8JBFnOwDYoeyIIELOb2b1LZPzmsxAMSmjZ7pR45bO2jh0O7O21O0WCqEcu76lt\/mVHUFlPQetzTiBIJ3rX\/9lz4pgPjAVfonzegRC00IIYlKxPY90fOwrYPeBC0ja70c9B63POIJhJpQ1aaJlsi9+yFAhuAMHFfg5anxWCfikUcD4XlUPb7Ab4FYIpwCVooaVQ5bQvUSC4AQQX+zlofVZLMG085aYAF5VD2+wG+BWCKTAkaKGlUOW0L1EguAEEF\/s5aH1WSzBtPOWmABeVQ9vsBvgVgikwJGihpVDltC9RILgBBBf7OWh9VkswbTzlpgAXlUPb7Ab4FYIpMCRooaVQ5bQvUSC4AQQX+zlofVZLMG085aYAF5VD2+wG+BWCKTAkaKGlUOW0L1EguAEEF\/s5aH1WSzBtPOWmABeVQ9vsBvgVgikwJGihpVDltC9RILgBBBf7OWh9VkswbTzlpgAXlUPb7Ab4FYIpMCRooaVQ5bQvUSC4AQQX+zlofVZLMG085aYAF5VD2+wG+BWCKTAkaKGlUOW0L1EguAEEF\/s5aH2OrCV4xhmN5L338qVNmzZpAyYKBbioHNpmN8CvEEyBQAgNCP7P\/zSSadNSKCCClygQ3ACCi\/2sEEwBSAjt7LNFLr+8jUIwBflF5RIXgeBimxWCKWikQlCtohSGTSQuUQjmfnorsnOCaglGQqfTqqSLQHCxzaGwBHntZa9evWTFihW+B23Lli1lypQp4n1vsO+L08xo5wRPP72RLFyYZmERudxF5dA2u2HxhwaC3pejl8YF3g187733yoQJEwKD4GmntZGzzhKFYGmdFeHfFYIKwVwM38i6wwrBXAyPYO+hEFQI5mIEKgRzIeUM3EOB4AYQXOznULjD6OiOHTtk8ODBMm\/ePKOyzzzzjPnbrVs387djx44yZswYqVq1agZUOr0iEJpagunJMApXuwgEF9scGgiOHTvW6MVtt90mtlJ9+\/Y1\/7eAbNCggfl\/0EkhqFZR0GMwW\/dXCAa0RIbosDcwYqHXtm1b6dKli+nvoIMh3kGnEFQIZgtCQZerEAwZBLt27Vq0N1chGKx6uKgc2mY3HnahcIcTWYIKwWDB5727AsENILjYzwrBFDij7rACIYVhE4lLFIIBusN+dowEuUtE5wR3yvfffy\/169eXKlUUgpEgWgqVVAgGBMEU+irQS9QSVAgGOgCzeHOFoELQ1\/BSCCoEfQ2UCGZSCAYEQb8HKKg7HJxWuagc2mY3HnahCIzEqjYLp5s2bVq0RpDfZ86cKYsXLw7FrhG1BN1QDoWgG\/0cOgjGLpexgMzlOkEgPHnyZHPr0aNHF4Mx3ykE3VAOhaAb\/Rw6CAIZCyH2D\/Mio9htdNl0Cr0WJztX4h3xpRB0QzkUgm70cyghaLfJ9ezZU7777jthz\/C0adOkWbNm2eSf2aM8atQo6dGjR4n3Ugi6oRwKQTf6ObQQzCrtEhSOKz58+HADwPvvv79Ud\/jII\/PljTcKpF69ekFUN6f3BAjr16+XunXrhuIkn1w0Xtsc\/IlN2e5nxvSyZctk4MCB5uSqIF6hG6rzBG2UGkFwWg3zkIMGDZJx48YVswytJUgH3X33cOnevXu2+yrw8nft2iUbN26UOnXqSOXKlQOvTy4qoG1Orp\/37RPhs3dvnhQU5MmePSK7duXJjh3lZOdO\/ubJ1q3lZNu2crJ9+y\/\/5nv72bMnT3bv\/uXatWsrmHIKC5Pr6XXrKvi+ID+\/eN7\/\/d+\/BQfBRMGQRK3JVJDEmsHch+U348ePN+7w0KFDDfTinWbjDYx469eoUYE0bFgghx32y7d0Xo0aIrt3i9gNFnv3ipQv\/8tvNg\/\/Lyj4b0n8xqdcOQbTL9\/bf8d+n5cnsmrVf8uzZXr\/xsow2UHlvR4oVKhQXsqX9z\/Q\/IzIZAaun\/L85IlVAD\/XaJ6yLYHAIehn25y3C7KxZtBCzx7ckAwEy\/bw0NapBIKVQNOmqd3fz3X5+fmSn\/+tLFyYF5wlmFrzsnMV0eE1a9YUHe46YMCA\/YIyXnc4O7UIttR4A6fg\/83VChUyawUm21I\/gzrZMr35sa5J3KewsMB4AzVr1pCKFfOMRV6x4i\/WPGLgL1Z+pUoiHHjOv\/lUq\/bLX77zfsjHTAIfyqEM\/lIuH+5t759OG9K51sVgkAZG4owY7zrBeJOlXgim416mM1hzfa2LyqFt1uhwLvQsVIERvw22EAxqDsFvPTOZT4HgBhBc7Ge1BFMghYXglVeKTJuWQgERvMRF5dA2uwF+hWAKQFIIuqEcCkE3+lkhqBD0JQEFghtAcLGfFYK+EFA8k1qCCoQUhk0kLlEIBnSeYCRGh6eSFoI9e4pMnRq12qdWXxeVQ9vsxsNOLcEUmGAhSFCE4IgLSYHgBhBc7OfQQrCk06azsVskGZApBBUIyYyXKOVVCIbcHWYR87x583JyrFZJA1chqBCMEtiSqatCMKQQZCvbkCFD4p7ynEwHZyqvQlAhmKmxFLZyFIIhgyCnxXCwaseOHc1e3rAkhaBCMCxjMdP1UAiGCII6J5jp4Z1eeS4qh7bZjYddaAMj6alsdq9WS9AN5VAIutHPoYBgsoeqZhdxpZeuEHRDORSCbvSzQrB05u2XQyHohnIoBN3oZ4WgQtCXBBQIbgDBxX4ODQT9HK8f9CJpSwu1BBUIvp4cEcykEAwoOqxzguHXFheVQ9vsxsMuNJZg\/\/79i97yFnYkqCXohnIoBN3oZ4VgCsRVCLqhHApBN\/o5FBBMgUOBXqIQdEM5FIJu9LNCMAWcKgTdUA6FoBv9rBBUCPqSgALBDSC42M8KQV8IKJ5JLUEFQgrDJhKXKAQDWiITidHhqaRCUCEYtTHrt74KQYWgr7GiEFQI+hooEcykEFQI+hq2CkGFoK+BEsFMCsGQQnDr1q1SqVIl8wlDUggqBMMwDrNRB4VgQBDcsWOHjBo1Sq644gpp3Lix1KhRw\/QvsClfvrz89a9\/lTPPPFPatMl9BeMNNIWgQjAbAApDmQrB3DMmr7CwsNBC8NRTT5X58+fLiBEjZM+ePXL\/\/febY\/WnTZsWSgjefbfIXXeFYehmvw4uKoe22Y2HXSiWyFgI9ujRQzZs2CAffPCBFBQUyHnnnScnnniijB8\/XiGYfc6VeAcFghtAcLGfQwXBTp06Sf369aVevXqyevVque+++6R3796hdYfVEgyYzFm+vYtAcLHNoYDgtm3bzJzgGWecIa+++qp07tzZzAfu3btX+vXrJ5MmTVJLMMsKX1rxLiqHttkN6zdwCG7fvr2QOcA333xTHnjgAalQoYIMGzZMRo4cKe+9954wEF955RU54IADzKdRo0ZmzvCggw4qTW+z9rsNjKglmDURh6JghaBCMBcDMW\/Xrl2FX375pUydOlXWr18vtWvXluOOO05WrVol3333nQwfPlxmzJgh7dq1k5NOOkny8vLkwAMPlHLlyuWifnHvoRB0QzkUgm70c+CWoDc6zBKZDz\/8UFgXiEt8++23y4ABA+Tll19Wdzgw5P9yYwWCG0BwsZ9DBUGiw1iC27dvlwULFsgJJ5wgLVq0kAkTJigEFYI5l4CLQHCxzaGCIJbgnDlzpHXr1lKzZk0THLn22mszBsGxY8fK5MmTjTKNHj1aunTpsp9iffXVV9KzZ0\/jipPivdxJ3WG1inJO5BzdUCEY8GJp4Me84FVXXWUiw0Cra9euMmvWrLQtwZkzZ8rixYtlzJgxwrrERO80AXDPPvusyVe1alWdE\/x\/CbioHNpmNx52obAECwsLZcuWLVKtWjUTHbZp06ZNxiJcuHChNG\/eXA4\/\/PCUnofexdjNmjUrsQxguWbNGrNTJVFSS9AN5VAIutHPoYCghc2+fftMUIS9w5mM\/vJKT6LMAJCteCW5w16XmXzPPPPMfnuWLQSvvjpfhg0rMIu7y3oCCFjpdevWTWghlzUZaJvje0JlqZ8Z08uWLZOBAwfG1fVctNXsHbY3srC666675OCDD87Y\/SmXl7tzAAMWHvN+gwYNknHjxhkw2oTFOHjwYGnbtq2ZLwR2RKfZu+zNZyFYq9ZEueWWn6V79+4Zq2tYC9q1a5ds3LhR6tSpI5UrVw5rNTNaL21z2e9nlt9Nnz7djJt4Bk9GB1SCwvJ27txZCIyWL18uu3fvlvz8fLMgOtGxWR06dDBzhl63ObZsa97yPYEN9h6zI2Xo0KEGZrGwS9TQRPlctASRBfu6sXqrVHHDTdI2l\/1+xhIkGDtx4sTgIJiMJcj2Oip7zTXXJGUpWpgRZMEaTBaC9joLS50TLPvKQV\/rnKAb\/RyqOUGvO8y2uI8\/\/tgsj+nWrZu0b9\/ewOvRRx81lmCy2+a8AY9Ebi7390aNyccc4ZQpU4pBVyHohnIoBN3o51BCkJ0iL774opmIZ52gDTwwfcgn1aCJN+hh\/f9Y8HnXCTZo0GC\/+UAsBIWgG8qhEHSjn0MFwZ9++kkefvhhueGGG0yEmHk\/9gqHLSkE3VAOhaAb\/RwaCHKE1pNPPlkq7zhmH+swyAilQtAN5VAIutHPgUNw3759hZwiza6QQw45RFasWCEcs\/\/WW2\/JhRdeKPfcc49ZrsIp0yxZ4LSZ66+\/PqnASKlkTTKDQtAN5VAIutHPgUNw69athawL5OxAIr9Lly41h6u+9tpr5iQZghIsomZLHctjXnrpJbnooouKXsaUJL8ykt1C8PHHf5C+fQ\/JSJlhL0SB4AYQXOznwCHIEhmCHQQq2I0Qawm+\/vrr5oh9osRYhqeffnrgvFBLUIEQ+CDMUgUUggEdoEB\/co7gokWLinUtS2I2b94sWIrsLSa6e+utt0rDhg2zNAT8FasQVAj6GynRy6UQDBCC8YYL+4gB4znnnGN2kKxcudLMBSoEc69cLiqHttmNh10o3GGr0lh+vGyddYAMQPsSdu\/vvI+Yd40EmdQSdEM5FIJu9HMoIMjGfIA3d+5cadKkiTRt2tTsDOnTp485Wp+9qmeddZZZKM2GZ17ElOisv1zAUSHohnIoBN3o51BAkPWBLIFZsmSJgSAvXuIQBY7Xf+6558wLlmbPnm2iwrx5TiGYC9QXv4cCwQ0guNjPoYAgW+SAH2+YA4DsC77gggvMkpmjjjpKbrzxRnPcDSfCcOIDx2CpJZhbELqoHNpmN8AfCghSiU8++cQc1cQSGc7+Y16QbXSc7cfe4bVr10r16tXNXt6bbrpJIZhbBuqJKo4cH+Yi+EMBwY8++kjeeOMNYeeIF4K87OjSSy8tAh6HHeA6KwRzTEA9VsqZMxQVggEtkeHkFj6sBbTu8LnnniuPPPKIeevbySefbMBH9BgI9uvXTy3BHHPQReXQNqs7nAs1M8frA0B2hnB0O3ODnCPIVrlDDz3UnCP4wgsvmKPdAeP8+fP\/c6T9LSWeLJ3timt02A3lUAi60c+hcIeZ+8MlBnR2icykSZOMxQcYOTjhwQcfNO5yq1atDBiDTApBN5RDIehGP4cCggCNNYC4vrVr1zYLpjlQlZ0h9jxB3GT2F\/OipCAjw9RVIeiGcigE3ejnUEAQALL+D1eXo7OIDPft21e2b99e9MIlXsJEpBgIlvSSpVxYiApBN5RDIehGP4cCghyewDtAWATNGsF4r97EFebUafYO886RIEGoEHRDORSCbvRzKCDInB+QYx0gif9\/\/vnncvTRRxc7QRo4fvPNN2YnSZDH7isE3VAOhaAb\/RwKCObChc3kPRSCbiiHQtCNflYIpkBHhaAbyqEQdKOfFYIKQV8SUCC4AQQX+1kh6AsBxTOpJahASGHYROIShWBA2+YiMTo8lVQIKgSjNmb91lchqBD0NVYUggpBXwMlgpkUggpBX8NWIagQ9DVQIphJIagQ9DVsFYIKQV8DJYKZFIIKQV\/DViGoEPQ1UCKYSSGoEPQ1bBWCCkFfAyWCmRSCCkFfw1YhqBD0NVAimEkhqBD0NWwVggpBXwMlgpkUggpBX8NWIagQ9DVQIphJIagQ9DVsFYIKQV8DJYKZFIIKQV\/DViGoEPQ1UCKYSSGoEPQ1bBWCCkFfAyWCmRSCjkOQt9717NnTvOvEm0aPHi1dunQp+kohqBCMIN98VVkh6DgEY0cJsBs7dqxMmTLFHOtvk0JQIeiLKBHMpBBUCBYNW170PnjwYOnatau0aVNcMApBhWAE+earygpBhWAxa+\/ZZ5+VMWPG7PeKT4WgQtAXUSKYSSGoEDTDtiQrkN8tBK++Ol+GDSuQevXqRXC4J1dllIN3QdetWzfw9z4nV\/PUc2ubq6YuvIhcyZhetmyZDBw40LzXPNbry0Uz8gp56XBAyR6rze1btmxZNPdHgOTee++VCRMmFJsLjJ0TrFVrotxyy8\/SvXv3gFqQu9uoh6EvAAAIZUlEQVTyBsCNGzdKnTp1ir0BMHc1yP2dtM2Vcy\/0HN9xxowZMn36dHNXJyGYSN68A3nNmjXmRe\/xkouWINbxhg0bjNVbpYob7rC2uez3M5bgnDlzZOLEiQpBL+yICDdt2rTYshjv7y7OCfJQ4InZo0cPIxsXkrbZjX7WFy3FaHNp84Gxc4K9euW7wANZt26dmTe5+eab5dRTT9U2l1EJuNzP6g4nMajz8\/OlceNGwpzgQQdNTOJKzaoSUAmEUQI82MePHy+NGjXKefUCDYyk01pAyEeTSkAlEH0JAL8gAIjkIgvB6He7tkAloBIIgwQUgmHoBa2DSkAlEJgEFIKBiV5vrBJQCYRBAgrBMPSC1kEloBIITAIKwcBErzdWCagEwiCByEHQe+Zg3759E+4qCYNwk6kDu2SGDBlSdIm3bSW1mYXlkydPlgYNGsi0adOkWbNmydw20LzxFsUnas+PP\/4ovXr1khUrVkjHjh2LHawRpTER2+bYMzS920ej3Ga73nfevHlmjMX2WZj6OVIQtIOC7XQMFo7aatu2bcKdJYFqeJI3T7RLpqQ2A87FixcbIACHeGcvJlmNnGW3SuA9MDdRe6pWrVqsr7mWxDiI0piI12Z2S8Q7LclCxI7vqLWZviRxGHJsW8LWz5GCYOzBCl5hoihRTSXtkimpzZMmTSraXuiFQRAncfiVvW0rlivJuz3S+yDwtqd58+bSv39\/GTp0qLF0vYft\/vDDD8UO2wjjmCipzYn2ydP+KLc5djx4+yXRuA2qnyMFwdinZqKTp\/0qZFjyed0e6uR1iRK1+eGHH5Zx48YVHTob+7QNS9tKqocXerEPAm97TjrppGKg48EwaNAg034g6LWkwj4mYi1+ax1aOdmtY7EPvyi3mbZZS\/amm24qdlhyGPo5UhCMfcqHfcD7BZF3gGPpMGB4zwpu7ty5c4tcXqxd22a2GLGXGJcQyy\/qEIy1ZL3tadKkSTFX3yuv5cuXx5VP7CsZ\/PZFtvPFA791eenbAQMGmLld4O6d3ohym716inyZ2403boPq50hBsKxagrGK5+epr5Zg2bAEvX0fBqso0w8BL9h5wIfR4o8UBMvqnGA8CNpDZUua84rinKC3rfFcQztHWJbmBEtqczwI8l6doObHMgnBRJ5a2OZ+IwXBKEUCkxlMsRau30hglKPDdp7IGxgJW9QwmT70mzcWAImCH1GPiJc0VRW2fo4UBBloUVoT5lcxyOddJ5jMOjhdJxitMVHSOsHYtZ5RXicYG\/BhjHvHta4TTIYOmlcloBJQCWRRApGzBLMoCy1aJaAScFACCkEHO12brBJQCfxXAgpBHQ0qAZWA0xJQCDrd\/Zlp\/O7du4UXpdesWVPy8vKSLrSgoED27t3rzPuUkxaQXpBVCSgEsypeNwr3LvGpWLGifPLJJ7Jly5ZijQeQxx57rFSoUKHY9wDwoYceMu+ZZodMlPeAu9HbZa+VCsGy16c5axHwe+edd+TMM8+UO++8U8477zw58sgjDQTZzrZp0yZjIfICHfb\/st939uzZxY4MK6my3j3UNp\/3iCaOG+vdu3exbVg2n98DJexSjdhlSTkTot4ocAkoBAPvguhWwAvBiRMnGiA1bNhQjjrqKMHCAzDHH3+8WR9mExbil19+KYceeqiUL1\/e\/Jt8QHT9+vVyxBFH\/Oc1qgfJ2rVr5bDDDpNatWpJuXLliq6P3XaVCHZ+IUjBiY6zim7PaM2TkYBCMBlpad4iCQC5t99+W95\/\/31p166dzJo1S4YNGybVqlWT119\/3Rx75U32kNjVq1fLHXfcYbaFXXrppTJy5Eg57rjjTNZFixYJ5wsuWbJE3nzzTeMeA0JvUgjqIMy0BBSCmZaoI+XFnohsm40Le+6550qdOnXMgbckDn5lzo+TQ0hbt26V1157TVq1aiWvvPKKLF261LxDGmuyRo0a8tZbb8nFF19sLMLYlAoEAa49ldpbnj22Si1BRwZtgmYqBN3u\/7RajzXInN\/8+fMN5Djai+jw888\/v1+5Xgjiqr744otywAEHmENVH3vsMbnuuuvk6aeflmOOOUa++OIL6devn3GZ\/UKQk7XjJQs672\/xjmSLd7pzWsLRiyMjAYVgZLoqfBVdtmyZcYOZA\/z+++9NYOSpp56SFi1aSPXq1YtVGAhyVt7UqVPliSeekB49epglNZyWU6lSJdm2bZvUr1\/fBFDat28v9913n+BCd+7cuVhEORVL0HvSdryN\/WoJhm9s5bJGCsFcSrsM3QsrEICdfPLJZjnM119\/bVrXuHHjovV+GzduNIES1hFaSxArkGU0AHDz5s3Gkvzss8+EvOeff74JiBx99NGyYcMGY1XaY\/it6NKBYOzhtbZMhWAZGpgpNEUhmILQ9BIxkdzp06dLnz59ZMGCBfLnP\/\/ZRHKx4AhwsPgZMHbo0EFWrVplIIiLy0nZ3377bTERfvrpp\/LNN98YSxJA2nTggQfKZZddZuYJ04VgSdFihaDbI1oh6Hb\/Z6T1QA1XlzeLEdDAWsM9xlJs3bq1ASFBjipVqsRdSP3uu+\/KypUrzRIbXGObKleubJbY8DcdCJb2ZkKFYEaGQWQLUQhGtuuCrziwI8rLCdcENi655BKzpm\/dunUyfPhwOeWUU6R27drGShwxYoR5gVS8lMwb4lJxh7lnt27d9ru1fd2nQjD4sRRkDRSCQUo\/wvdmmQvr\/Xbt2mXeHkYkd9++fSYy\/MgjjxjXt1OnTiao8dFHH8lf\/vIXkz\/espd0IJgJESoEMyHF6JahEIxu3wVec4IjsXuBASEBjdiDFOLlTaUBJb2jOZXyuEYhmKrkysZ1CsGy0Y\/OtCJ277BdgJ2qAHTvcKqSKzvX\/R8GDkwTJt51hQAAAABJRU5ErkJggg==","height":194,"width":321}}
%---
%[output:4faf0d88]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAUEAAADCCAYAAADXXXDzAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQ1wVsXVPlgqyWfFEcoQEAVpoT\/MlFqsBtpSqy1TqOhYf\/iZKZhSBBXMDBAIFCiSCokBxwSUBg0ROg2ilTGVWqm\/tQ7TjKUtU0WQT8nIjyAD1lJLqPTjm3Mz+7LZ3PvevXf37u6977kzGUjes7vnPGfP8579vd3Onj17FughBAgBQqBAEehGJFigniezCQFCwEOASJA6AiFACBQ0AkSCBe1+OeP\/9Kc\/weTJk6G5uRlKS0vlCqVMat++fVBWVgazZ8+GCRMmpEx7UlcFASJBFfQKpCyRYIE4ukDNJBIsUMdHMZtIMApaJJs2BIgE0+YxC\/oSCVoAnZo0hgCRoDGo09uQHwm2trbC3Llz4aqrroJly5ZBz549fQ08deoUbN68GX71q1\/B\/v37oXfv3vDjH\/8YpkyZAhdccEGuDO7UevPNN+Hhhx+GF198EfD3ESNGwF133QXf+MY34LzzzvNkP\/74Y3j88ce9+UlW38SJE736+vTp48mcOHECpk2bBrfccov3\/7Vr18Kll14Ka9asgS996Utw7NgxWLduHTz11FPQ3t4OY8aMgZtuugmWLl1Kc4Lp7aaxNScSjA1d4RQUSfCvf\/0rzJs3D4YPH56XAP\/1r3\/Bz372M3jmmWfghhtugHHjxsHu3bth48aNMHr0aLj33nvhM5\/5jAfka6+9BnfffTcMHjwYpk6dCkVFRfDkk0\/CH\/\/4R6iuroYf\/vCH8MEHH3jtIgHjIsbIkSPh73\/\/u1ffgAED4IEHHoDLL788R4L\/\/Oc\/vd+RJA8fPgzjx4+HTz75BO655x545513vHa+\/OUvw7PPPgvbt2+H06dPQ1VVFS2MFE7X9iwlEiwwh8cxlyfBiy66CGbNmgVf\/epXO5GYX71IfuXl5bBy5Uq47bbboFu3bp0IDzNJzOA+\/PBDjwCREFevXg0XXnihJ4ckumTJEvjUpz7lZWlPP\/001NbWepndt7\/97VyTe\/fu9TK\/sWPHwoIFCwDJD3\/HfxsbG2HQoEE52U2bNnltPPTQQ\/DNb37T+\/uZM2e8LBF\/UFdaHY7TS9Jbhkgwvb4zpjkjwcWLF3vDWsyekCxYFuenCGZVixYtgrfeessjon79+uXEcIhcWVkJ+C9mb7g9BckQM74f\/OAHvnYhIc6ZM8f7DMvwbSOJYQb3+uuve2316NHDI8HLLrvMq7O4uDg3lMZ28cG\/88Pxt99+28suMUskEjTWtZxoiEjQCTe4rQQjQaYlDjHXr18Pn\/vc5wIVx+zujjvugO7du8P06dPh\/PPP7yS7ZcsWeO+99+DRRx\/1hqZh+xDZPB\/uU8RsT3ywPszkmpqavHlHJMGvf\/3rsHDhwlwGmq8O9hkOnYkE3e6PurUjEtSNaAbrYySIQ1okCJyXu+KKK7zsi2VZotmMVHbt2hWISP\/+\/T3SOn78uDIJ4mIJDpN5EhQJk0gwg51Tg0lEghpAzHoV\/Jzg1VdfDU888YQ3V8cWLNB+lGloaPCgwOEt\/n\/Dhg3wn\/\/8p8vQU8QLF1qChsNYxwsvvOAtwKxatSrScFgkQXEYzg+p6cRI1ntxsH1EgoXre2nLxdXhkydPettj2travGExLjz4baPxW4TARnGVF+fecEvLfffdB\/\/973+9hREcxvJzjUeOHPG2rODwG7NOHPLmWxi59tprAectcf4Qh8N+Q2fUSawDt+Pg33G1mhZGpLtFZgSJBDPjyuQM8SM4zN5uv\/12uP76672VW\/y9pqbGW5jo1auXpwzbIrNt27bclpb333\/fI7N333230wqtuEUGy+PWl4MHD3pEi8PvfFtkLr74Ym84PHTo0NwWGT8SZDphdolzlbRFJrl+k5aaiQTT4imLevqRIK7I4jYT\/HnwwQc94hNJEFVmm5txY\/KePXu87S\/f+ta3YObMmTBs2LDcooW4WRrLXnfddV7GyS\/ARNksHbSIgnVg5odD7Y8++sjLGHFluL6+3ttTSAsjFjubhaaJBC2AnsUmkSj9SDCLtpJN2UKASDBb\/rRmDZGgNeipYUUEiAQVAaTiHQgQCVJPSCsCRIJp9RzpTQgQAloQIBLUAiNVQggQAmlFgEgwrZ4jvQkBQkALAkSCWmCkSggBQiCtCBAJptVzpDchQAhoQYBIUAuMVAkhQAikFQEiwbR6jvQmBAgBLQgYJ0E8VYAH7tnRJHZ7B15\/jg9e2c6fP9ViJVVCCBAChEAAAkZJEAkQr1jib+rATbb4Ih7+BmA\/XfEgPf7QQwgQAtlDAN8Rgz82HiMkyO5xw0s08eEzQbxRBK9k8rstmAGC5FdRUeG9YIceQoAQSCcCZ84MgH\/84244eXKCrwEHDhy0QoRGSJC3WBwOs+yQyeCrFPFWD\/5ht5jcfPPNcOWVV0Lfvn3T2Qt8tN65c6d3E8vy5cutdICkgST7kkY42fpV\/VdX1wO2bbsUkAD9nu7dO0Z3\/fpNhieeuL9L7CdrXUftVkmQZYijRo3y5giR7PDqdrwifciQITn7+Xdc4GsS8RbirDyHDh2CrVu3eq+UvOSSS7JiVs4Osi\/dLo3qv927z4eqqt7Q2loUaPiAAWdg27ZD0LPn\/3lXmuG9kfj4JUAm0LNKgqKBIimyzxkJ4o3A+ELukpISE9gYaQNtPnr0qGcTvms3aw\/Zl26PyvjvmmsgL+nhG09ras7ADTec6QIG3h7e0tICdXV1RIKIDiPBSZMmdUqL\/S71THfXOqd9e3s74G3L+ErKLJIg2Zfunurnv5\/8BKCxMdguJL3x4wHq6+Vstx3fVjNBfPsX3hyM76fF4W\/QdUy2QZJzZTwpIol4uLlSKuv+27OnHV599SP45S97w2uvdfeFHUkPf5qaOv6N+tiOb6skiGDx+wTZKxj5+UCUsQ1SVKdGkc96EJF9UXqDO7Lf+Q7AK6\/kz\/ZefBFg8GB1nW3Ht3ESjAOZbZDi6CxbhkhCFik35bLgv7Y2nLMD+MUv8pNeRQXAXXfp94Pt+CYS1O\/TSDVmIYjyGUz2ReoOxoT37gWYOTM422ND3N\/9Lvk5ayJBCbfbBklCxdgiRBKxoXOiYJr8JzPEff55gM9\/\/hy0JuyzHd+UCVoOJROdzKaJZJ899N94A2D27PzZ3s03A6xaFayjCf8RCUr0EdsgSagYW8REJ4utnIaCZJ8GECNUkS\/bYyu3L78sv4prwn+245sywQgdLAlRE50sCb1l6yT7ZJGKJ3fgAAAeoApayUXiw2xwzpx49ZvwH5GghG9sgyShYmwRE50stnIaCpJ9GkDkqsCV3KoqgA0b\/OuNk+3l09CE\/2zHN2WCevto5NpMdLLISmksQPapg4nEV1aWP9ubMQNg4kT5Ya6sVib8RyQo4Q3bIEmoGFvERCeLrZyGgmRfPBBrawHmzw8uixlflLm9eFoAmPCf7fimTDBu79BUzkQn06RqrGrIPjnYZLK9738fYMEC\/dkeDYfPnj0r5yZ7Ura\/KZK0nEgiSXSTr1vFfzLEt20bwLBhydsR1IKKfbJa245vygRlPZWQnIlOlpDqUtWSfedgQtLD87Z4C4vfo3oRgZRDIgqZ8B+RoIRTbIMkoWJsEROdLLZyGgoWun1IfMuXd9ywEkR8uGF51iyzw1xZ15rwn+34pkxQtjckJGeikyWkulS1hWifzDC3pQXgK1+RgtCqkAn\/EQlKuNg2SBIqxhYx0cliK6ehYKHYV1LSD8aNK8q7adnEaq4Gl3WqwoT\/bMc3ZYK6e03E+kx0sogqaRXPsn2Y8eExNfxXfFyc34vjWBP+IxKU8IxtkCRUjC1iopPFVk5DwazZh+8Euv324Pm9W2\/tuHMvzg3LGuDWXoUJ\/9mOb8oEtXebaBWa6GTRNNIrnXb7MMtbvx5g5Up\/XPDNaUuXHocf\/egiekdMzK5DJCgBnG2QJFSMLZJ2kggzPI32IfHdfz\/AunXBGd\/27QBDh5o5URGGcZKfm\/Cf7fimTDDJHiRRt4lOJqFGYiJpsQ+J7+c\/D36LWtAxtbTYF9fBJuwjEpTwjm2QJFSMLWKik8VWTkNBl+2LS3w8LC7bp8F9dHZYB4g66iAS1IGinTpcIwmZq6iibGVxzT7dXjZhn+34puGw7l4TsT4TnSyiSlrFXbAPiW\/FCoBHHgme44tCfJQJau0i1l+pSySo15+Ra3OBJCIrHaGALfuQ+OrqAB58MJj49u+PYEiAqC371DWXq8GEfZQJSvjCNkgSKsYWMdHJYiunoaBp+3BFN+jduEncwWfaPg0uiVSFCftsxzdlgpG6hH5hE51Mv9byNZqwr7UVoLRU\/1BXxkoT9snokZSMCfuIBCW8ZxskCRVji5joZLGV01AwKfvwFszBg4OPrMWd44tqclL2RdUjKXkT9tmOb8oEk+o9kvWa6GSSqiQiptO+fLez4FD3uecAvvCFRMwIrFSnfWY1l2vNhH1EghK+sA2ShIqxRUx0stjKaSioal8Y8dXUANx2mwZFY1ahal\/MZo0VM2Gf7fimTNBYd\/JvyEQns2liXPvuuw9g8eKummPGd\/31AHPnunFJQVz7bPokStsm7CMSlPCIbZAkVIwtYqKTxVZOQ8Eo9oVdTWVqni+K2VHsi1KvK7Im7LMd35QJWu5tJjqZTRPD7Asb7rpIfDyeYfbZxF5H2ybsIxKU8JRtkCRUjC1iopPFVk5DwSD7Zs8GWLvWf7j75JMAV16poXEDVRSq\/3RCazu+jWeCNTU1MGjQIJgwYUIOR\/xbQ0MD9O\/fH5qammDIkCGdMLYNkk6Hi3UVUhAdOVLkexMzzvN973sAixa5Mc8Xxd+F5L+ioqIo0EjL2o5voyTIyG7lypU5EtyyZQvs2LEDqqurYdeuXYAyjY2N0KtXrxyItkGS9mYMwawH0Z497d5NzK2tXQMoiRMcMVygVCTr\/jNhn+34NkKCp06dgsrKSi\/Tw4fPBPnM8MSJEzBt2jRYsGABlHJHAGyDpBQlIYVNdLIk9Q+q+w9\/ALjmGv\/h7tatAFdcYUMr\/W1m1X8MKRP22Y5vIyTIdz2e9Bg5Tpo0ySM99vuoUaM6DZcZSOXl5XDjjTdCSUmJ\/t5sqUbsZEeOHIG+fftCcXGxJS30NIuLHHfe2fWNayUl7TBwIMBjj6VvuBuGTJb852dr0vZh39+5cydUVFRAc3Nzp+QnDHtdn1slQTHzCyNBNHrq1KkwZcoUXfZbr+f06dNw7Ngx6NOnD\/To0cO6PnEU+P3v\/wdmzuzbpSi+f+P55\/839fblwyQL\/rNp36ZNm2Ajvr0KoDBJMGomWFtbCyNGjMhUJogYHD161LMpqYnnOMQmUwaHu+JcH87zjRsHsHp1u1dFmu2TwYDsk0EpWAYzwZaWFqirqytMEkRoaE6wHd5\/\/33o169fKkhw716AL37Rf67Pb0+fiTkltTBUK032qeGHpQt6ThABoNXhdJAgvmT8lVc6d3jM+saOBXj44eBAIJJQJwmbNZjwX8GTIMsGaZ+gm5ng5Zd3va4KyU\/2VmYTQZR1ksi6fQVHgnEcahukODrLlnGRJJYs6Xj9pPjgHCAOeaM8LtoXRf8wWbIvDKHwz23Ht\/HV4XBIukrYBimOzrJlXAmioDO8mPX9+c8AvXvLWtRZzhX74mkfXorsC8coTMJ2fBMJhnko4c9dCCLVIW8+iFywL0kXkn3q6BIJSmBoGyQJFWOL2AqibdsAxo\/vrDZmfWPGADQ0xDanS0Fb9umzIH9NZJ860rbjmzJBdR8q1WA6iMrKOk5u8A+S365dAD17KpniW9i0ffotIBJMegsXkaBEr7UNkoSKsUVMkUTQFpek7+szZV9sBygWJPsUASzEfYJxICMSjINaR5kk5\/tktCKSkEHJXRkT\/rMd3zQcttz\/kupk3bp1HfLecgtAba1Zg5Oyz6wVwa2RfeqeIBKUwNA2SBIqxhbRGUSHDwNccklX8nvppY6M0Maj0z4b+oe1SfaFIRT+ue34pkww3EeJSugIot\/+tuMNbPwT5VRHkgbqsC9J\/VTrJvtUESzAs8NxILP9TRFHZ9kyKkH00EMAs2Z1Jb+kFztkbUM5FfuitGNLluxTR952fIdmgnjn39\/+9je49tprfa3Fq4Tq6+th+vTpna7EV4fmXA22QdJpi1hXnCBauBCgutrNzE+HfUnirbvuOP7TrUOS9Zmwz3Z8S5EgXnk\/ceLETrc9I\/BM+eHDh3d5L4hOx9gGSactKiRRUwNQWZkO8mNamgiiJP0TVjfZF4ZQ+Oe24zuUBHmyYy9IYjdC44uRZsyY4b0TJMnHNkhJ2iYTRPX1AOXl6SI\/IsEke425umX6p6o2tuNbigR5ImQGJ5398cDaBknVyfnK5+tkv\/41wK23ppP8iAST7DXm6iYSFLBmZMS\/MtOEOwqNBA8cALjssnSTH5GgichIvg0iQR+MbRCSjTaT714dLYidzG+Ts+wFpqZ0jtKOiSCKoo9uWbJPHVHb8e07HObn\/GRMTHpobBskGQziyrAgGj36Ujh4sHuuGlf2+cW1izJBVeTcKG+C5G3Ht\/ScoE2X2AYpSdvFs71ZIT8iwSR7jbm6iQTNYZ23pSySoHirC5Lfvn0A3c8lg46gr6aGiSBS01CtNNmnhh+Wth3flAmq+zBSDRs2AEybdq4IvqC8ouIE3HFHz1S8cjOSsXRiJCpczsmbIHkiQQm32wZJQkUpEb9Fj7feSscrN6UM9BEyEURxddNRjuxTR9F2fFMmqO7D0BryzftREIXC57QA+U\/dPUSCEhjaBklCRV+Re+8FWLbs3Ed+ix4URHHRdaMc+U\/dD7bjmzJBdR\/61iAOfZ9\/HuC73+0qSkGUkAMMVUv+UweaSFACQ9sgSaiYE4m65YWCKAq67smS\/9R9Yju+KRNU96FXw3PPAYwdm3\/o69cUBZEmB1iqhvynDjyRoASGtkEKU1HM\/qqqABYvDivV8TkFkRxOrkqR\/9Q9Yzu+KRNU8CG+tGj+\/OjZH98kBZGCAxwoSv5TdwKRoASGtkHyU5HP\/nDV9\/XXAT77WQljBBEKouiYuVSC\/KfuDdvxTZlgRB+uXw8wY4Za9keZYETQHRYnElR3DpGgBIa2QWIqitnfG28AXHCBhAF5RCiI1PCzXZr8p+4B2\/FNmaCED8+cAfj0p\/Vlf5QJSoCeEhEiQXVHFTwJ7tu3D8rKyuAwvjkcAPzuJrQJkrjyu2MHwMiR6o5nNVAQ6cPSRk3kP3XUbcY3am89E0QANm\/eDNXV1VBcXOyLqC2Q+FMfSd3zR0GkHkQ2ayD\/qaNvK76Z5tZJcMuWLdDW1pb3jXWmQfrwQ4BevZIZ\/opdhoJIPYhs1kD+U0ffdHyLGlsnwZqaGmhoaMjp1dzcDKWlpZ30ZCCVl5fDjTfeCCUlJerIB9Rw551F8Nhj5z78+ON2OO+8xJrzNksfOXIE+vbtG5gJJ9d68jWTfcljnGQLSfsP+\/7OnTuhoqIC\/GI\/SducyARPnToFlZWVMGrUKO\/F7kh28+bNg6amJhgyZEjOfkaC+IepU6fClClTEsGGf88HXnb66qsHEmmHr\/T06dNw7Ngx6NOnD\/To0SPx9kw3QPaZRlxve0n7b9OmTbBx40ZP6YIkQdFdIimyzxkJ1tbWwogRIxLJBIuLi3Lq4PwfXnZq4kGbjx496tlUVHROBxNtm2iD7DOBcnJtJO0\/zARbWlqgrq6OSBDdyEhw0qRJnYbESc4ZfPIJwPnnn+tEL70EgO\/\/MPXQnJIppJNph\/ynjmuS8S2jndU5QXy159y5c2HRokXe8BfBwDnCxsZG6MWtTCQFEn\/2N6nV3zAnUBCFIeT25+Q\/df8kFd+ymlklQVSS3yfYv3\/\/LvOBKJMESPzb3mwRINpGQSTbVd2UI\/+p+yWJ+I6ilXUSlFFWN0ji8bf9+2W0SEaGgigZXE3VSv5TR1p3fEfVqOBI0CUCpEwwand1T55IUN0nRIISGOoAqa2tY8ED\/8XH5hCYN5mCSKIDOCxC\/lN3jo74VtGiIDJBkQBnzQJYs0YFNn1lKYj0YWmjJvKfOupEghIYqoAkEiC+BnPpUolGDYlQEBkCOqFmyH\/qwKrEt3rrDlygIGOECkj8HODq1QBz5si0aE6Ggsgc1km0RP5TR1UlvtVbzzgJ8gS4ZAnA8uU6INNbBwWRXjxN10b+U0ecSFACwzgg8fsAZ84EWLdOoiELIhREFkDX2CT5Tx3MOPGt3uq5GjK5MHLddQB4\/A0fV1aBg5xGQaSzO5uvi\/ynjjmRoASGUUDavRtg2LB0ECBqSUEk0QEcFiH\/qTsnSnyrt9a1hsxlguw2aNczQOYKCqIkurW5Osl\/6lgTCUpgKAtS2giQMkEJ5zsuQiSo7iDZ+FZvyb+GzGSC\/Erw2bNJwaW\/Xgoi\/ZiarJH8p442kaAEhmEgLVsGgJug8cGXo0+fLlGpIyIURI44IqYa5L+YwHHFwuJbvYX8NWQiE0zjMJjmBJPu2mbqJxJUx5lIUALDfCClmQBpTlDC+Y6LEAmqO4hIUALDIJB++lOAFSs6KkjTPCBvMgWRRAdwWIT8p+4cIkEJDINASnsWSJmghPMdFyESVHcQkaAEhn4g8cfi0poFEglKON9xESJBdQcRCUpg6AcSywKrqgAWL5aoxFERCiJHHSOpFvlPEqg8YkSCEhiKILE9gWk5FZLPRAoiiQ7gsAj5T905RIISGIogsSwQL0wdOFCiAodFKIgcdo6EauQ\/CZBCRIgEJTDkQVq4sBReecX922EkzPJEKIhkkXJTjvyn7hciQQkMeZBGjiz1SmzdCnDTTRKFHRehIHLcQSHqkf\/U\/UckKIEhA+no0Wb49787SDDNK8K8yRREEh3AYRHyn7pziAQlMGQg7d\/\/rid9zTUAL78sUTAFIhREKXBSHhXJf+r+IxKUwFAkwaxkgTQnKOF8x0WIBNUdRCQogSGChJuj29tLnb8uX8KcTiIURFERc0ue\/KfuDyJBCQwRJLYgkqWhMGWCEs53XIRIUN1BRIISGPIkmKWhMJGghPMdFyESVHcQkaAEhkSCEiA5KkIk4ahjJNUy4T8iQQlnMBLMwjE50VwTnUwC4sREyL7EoDVSsQn\/EQlKuJJIUAIkR0VMBJFN08k+dfSJBAGgpqYGGhoaoH\/\/\/tDU1ARDhgzphCwjwTFjALZvVwfdpRooiFzyRnRdyH\/RMRNLFDwJbtmyBXbs2AHV1dWwa9cujxAbGxuhV69eOawQpLFjW6Gl5WYYPXqAOuoO1UBB5JAzYqhC\/osBmlCk4EkQSW\/QoEEwYcIEOHHiBEybNg0WLFgApaUdx+PwsQ2SupuDa2hra4ONGzfC1KlTPRyy9pB96faoCf\/Zjm+rb5s7deoUVFZWwqRJkzzSY7+PGjXKI0WRBMvLy+Hqq69Od68StD906BBUVFRAFm1DU8m+dHdXE\/5jbTQ3N3dKfkwhZ5UExcwviAQPHjzoEUVra6spXKgdQoAQMIgAJje1tbUwYID56S6rJCibCaIvkAjxhx5CgBDIHgJIfjYIEJG0SoKogMycYPZcThYRAoSAKwhYJ0GZ1WFXwCI9CAFCIHsIWCdBlg3m2yeYPdjJIkKAEHAFASdI0BUwSA9CgBAoPASIBAvP52QxIUAIcAg4T4L79u2DsrIyOHz4MMyYMcPbSJ22B+c9Fy5cmFObtyOffWHHCV3AgV\/YYvoE6c22ROHJoPHjx3unhIqLi71irvpZtI\/XE\/UePnx47oRTWuxjuzKeeeYZD3vRF1nyn0yMOE2C\/D5C7Gy4sVrcSC1jpG0ZP6JAnfLZl4YFIxYsK1euzG1uD9IbyY73H5bFB7\/UXPWzn314umHz5s2dCBztEPe4umwf+ggfPJAg6p0l\/8nGvdMkiN+6K1asgNWrV3tniXkHsQxC1lBbcuJeSF6PfPbV19eHHie0bRNeeIEPO\/aI\/w\/a8jR06FCYO3cuLFq0yLsgA8mEnRM\/fvy4U35mPvOzD\/sgHiUTRyRI5GmxT+w3fFwF9bs0+S9qXDhNguK3Lh84\/AULUY02Kc8PkcThU5B9a9euhfvvvz\/0OKFJO4La4kkv3+b3r33ta52IDr8A5s+f79mJJMhnVy75WcziWXbI8GBHvcQvtLTYx7648N977rkn8BhrWv0nEyNOk6CY+bkUHDLgogwfDJgBYRDh\/CbOh\/3mN7\/J3aCDmS2zD48P4TFBdpFE0HFCWR2SlMu32Z3Xe+DAgZ1uCOJx+ctf\/uKLg3ibUJJ2RCF5NiWD\/po3b553\/RsSOX8DUlrs42MKMeAvMMmC\/2T6jNMkmIVMUHSCTIZAmeC5obLtjD9oPlecB0xjpsSTOH5BZzGTTz0JZmFO0I8E2Txnvrkwl+cEeZv8hot+V6OldU5JhgTxFqS02Rc0qsrCnK4M8fEyTmeCrq4aRgFZzGZlVw3TsDosLoTg71lbXRRJIWjxI02r3\/mmlbLmP5lYdZoE2ZxalvYJRtkfR\/sEZbpwsjL59gmKr4NIyz5BcXEHEeT7Je0TTLZPUe2EACFACDiFgPOZoFNokTKEACGQOQSIBDPnUjKIECAEoiBAJBgFLZIlBAiBzCFAJJg5l5JBhAAhEAUBIsEoaDkq67fax1TlbzlxVH1n1OLPBYsbiXkl+bPFMrcaieeKnTGYFPEQIBLMQEfgj+Kl5WIJ12AXN+brJEG0Nej2GddwKER9iAQz4HUiQTUn+p3N1k2C+W4TUtOeSqsiQCSoiqAD5WVIkA3JRo8eDVVVVcBv9BUvffV7CTY\/5MZLYfFo3I4dO7yLIDDA+YP3CEnQkJGvR9xszE4yTJ8+HWbNmpVDVtQn6GYeLIB6TJw4MXe\/IcvC2EUHeEZWfMRLLsLKiLaJ+LH6RfvSeBWcA907cRWIBBOHOPkGZEkQCQIf\/nYWDMw1a9Z4N6EgQbCbk2fPnp0jEqwfCYqVY0TGThnIkqBYD9Y5efJkYCTHfhdPL\/Bt++mHNjz++OOefo888kjulh7kmV3FAAACl0lEQVQ2NRCGjx85qWSCjKRLS0s73TvoR7bJ9w5qIQwBIsEwhFLweViQowksMPksye9vKMuTCl7ygMcWV61aBRjUfJaH\/5fNBBl58fVgeV53vHZfzNjEcmHZlEheQTYytwZlrIyQ87nf73UPrD6GDT9Hy5+FZ1imoHtlXkUiwQy4OGh1mB+O+QVgULYjXvfF35PH4OLJSyYT5ImVvx6LP8z\/9ttvd7qTD9viSZC9YgHtClqVFUkvX0YX9OUQdTjMdyExs+Y\/c\/leyAyEQWwTiARjQ+dOwSiZILuolQU6Dkf9HkageOEpG2ry5MVnZLIkyL9sim+TbeORJcGw98yI18WzS2z9Vs6DMsU4w2FxeC\/iGnVrjTs9LNuaEAlmwL8qJJhvwYARZZKZIA+\/3xVPUTNBlj3i1f1I+Ki7uFDCt6mLBIPmASkTdD\/AiATd91GohnFJMGieTnwJkuycoN98I1scCMqswu5NFHWUsZVlXBdeeCG8+eab3ntM\/FaFEdh8c4JBXxB+ZWT0ojnB0K5sRYBI0ArsehtVCUBxDssvMxLn89iWELaKi9bg6zTxYe8SZvOU\/OKBuDosElxYJoiE6rc67Eewoo75NpGrrg7nmwfkPU2rw3r7va7aiAR1IWmxHhUSRLXFfW78e4SZWbwMkh9mWSdPnsyRnrh3D+vAV1Piwy9iiIs4\/B5AGRLE+sS2xP14bEiMGSy\/1SfIRSr7BHFPI249wpVtv4e3L2xl22IXKuimiQQL2v3xjZch3vi1q5cMWxU2PVdHJ0bUfZpUDUSCSSGb8XpdJ0HUT8xC87lEPDus2310dlg3ovrqIxLUh2VB1eQqCUaZCxQdxt8io9OZdIuMTjT110UkqB9TqpEQIARShMD\/A6hi8tffTTmiAAAAAElFTkSuQmCC","height":194,"width":321}}
%---
