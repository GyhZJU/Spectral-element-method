% 目的：Flexural wave propagation in beams with periodically attached vibration absorbers: Band-gap behavior and band formation mechanisms
% Figure 3
% 
% 状态：成功复现FRF曲线
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
geometry.ucL = 0.1;      % 单元长度 [m]
geometry.N = 6;          % 单元个数
geometry.L = geometry.ucL * geometry.N;        % 梁总长度 [m]


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
omega_values = linspace(1, 12000, 1000);
fre_values = omega_values ./ (2*pi());  % 转换为Hz
T_N_results = zeros(size(omega_values)); % 存储透射率结果

%% 主循环
for i = 1:length(omega_values) %[output:group:9c920127]
    omega = omega_values(i);
    
    % 计算当前频率下的透射率
    T_N = calculate_transmission_rate(omega, material, geometry); %[output:1c1d99d2] %[output:84b53c01] %[output:4171ae8c] %[output:5651274e] %[output:2e8e965b] %[output:5d4f7ab4] %[output:878ee9df] %[output:7aa25792] %[output:8cfd2c89] %[output:5b85e6c2] %[output:2c863765] %[output:1dfd022c] %[output:75481f82] %[output:5664c51f]
    % 转换为分贝(dB)：T_N(dB) = 20*log10(|T_N|)
    % T_N_results(i) = T_N;
    T_N_results(i) = 20 * log10(abs(T_N));

end %[output:group:9c920127]

%% 绘制结果
figure; %[output:3c2354af]
plot(fre_values, T_N_results, 'b-', 'LineWidth', 2); %[output:3c2354af]
xlabel('频率 [Hz]'); %[output:3c2354af]
ylabel('透射率 T_N [dB]'); %[output:3c2354af]
title('透射率随频率变化曲线'); %[output:3c2354af]
grid on; %[output:3c2354af]
xlim([0 2000]); %[output:3c2354af]



%%
function T_N = calculate_transmission_rate(omega, material, geometry)
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
    geometry.ucL;
    % La和Lb段单元长度
    La = geometry.ucL/2;  % 假设La和Lb各占一半长度
    Lb = geometry.ucL/2;

    % 计算单元矩阵（La和Lb段）
    D_La = spectral_element_matrix(omega, material, geometry, La);
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

    % 构建广义特征值问题矩阵
    A = [D_RL_, D_RR_;
         zeros(2), eye(2)];
    B = [-D_LL_, -D_LR_;
         eye(2), zeros(2)];

    % 求解广义特征值问题
    [eigenvectors, eigenvalues] = eig(A, B, 'vector');
    eigenvectors = eigenvectors (1:2,:);
    mu = log(eigenvalues);

    % 识别互为倒数的特征值对
    pairs = zeros(4,1); % 存储配对索引
    tolerance = 5e-3;   % 容差，用于判断是否为倒数
    

% % 寻找共轭对
% conjugate_pairs = {};
% remaining_values = mu;
% 
% while ~isempty(remaining_values)
%     current_val = remaining_values(1);
% 
%     % 寻找其共轭复数
%     conjugate_val = conj(current_val);
% 
%     % 在剩余值中查找共轭复数
%     idx = find(abs(remaining_values - conjugate_val) < 1e-10);
% 
%     if ~isempty(idx)
%         conjugate_pairs{end+1} = [current_val, conjugate_val];
%         remaining_values([1, idx(1)]) = [];
%     else
%         conjugate_pairs{end+1} = current_val;
%         remaining_values(1) = [];
%     end
% end
% 
% % 显示分组结果
% for i = 1:length(conjugate_pairs)
%     fprintf('第%d组: ', i);
%     disp(conjugate_pairs{i});
% end

    % % 寻找成对特征值
    % for i = 1:4
    %     if pairs(i) == 0 % 如果该特征值还未配对
    %         for j = i+1:4
    %             a=real(mu(i)); b=real(mu(j));
    %             ab = abs (a + b);
    %             c=imag(mu(i)); d=imag(mu(j));
    %             cd1 = abs (c + d);
    %             cd2 = abs (c - d);
    %             if ab < tolerance && (cd1 < tolerance || cd2 < tolerance)
    %                 pairs(i) = j;
    %                 pairs(j) = i;
    %                 break;
    %             end
    %         end
    %     end
    % end

    % 寻找成对特征值
    for i = 1:4
        if pairs(i) == 0 % 如果该特征值还未配对
            for j = i+1:4
                current_val = mu(i);
                conjugate_val = mu(j);
                result = current_val + conjugate_val;
                if abs(conjugate_val + current_val) < tolerance ...
                        || (imag(conjugate_val + current_val) == 2*pi())
                    pairs(i) = j;
                    pairs(j) = i;                
                    break;
                end
            end
        end
    end

    % 确保找到两对互为倒数的特征值
    if sum(pairs == 0) > 0
        fprintf('\n=== 特征值配对错误 ===\n');
        fprintf('当前频率: omega = %.3f rad/s (%.3f Hz)\n', omega, omega/(2*pi));
    end
    
    % 重新排序特征值和特征向量
    sorted_indices = zeros(4,1);
    count = 1;
    % 找出能量流大于0的特征值（正波）放在前两位
    for i = 1:4
        phi_u = eigenvectors(:, i);
        phi_f = (D_LL_ + exp(-mu(i)) * D_LR_) * phi_u;
        phi_f_conj = conj(phi_f);
        energy_flow = 1i * omega * (phi_f_conj' * phi_u);
        if (isreal(mu(i)) && real(mu(i)) > 0 && count <= 2) ...
                || (real(energy_flow) >0 && count <= 2) ...
                || (imag(mu(i)) == pi() && real(mu(i)) > 0 && count <= 2)
            sorted_indices(count) = i;
            sorted_indices(count+2) = pairs(i);
            count = count + 1;
        end
    end

    % 重新排列特征值和特征向量
    eigenvalues_sorted = eigenvalues(sorted_indices);
    eigenvectors_sorted = eigenvectors(:, sorted_indices);
    mu_sorted = log(eigenvalues_sorted);

    % 分离正向波和反向波特征向量(位移)
    phi_u_plus  = eigenvectors_sorted(:, 1:2);   % positive wave
    phi_u_minus = eigenvectors_sorted(:, 3:4);   % negative wave


    % 计算特征向量(力)
    for j = 1:length(mu_sorted)/2
        phi_f_plus(:,j)  = (D_LL_ + exp(-mu_sorted(j)) * D_LR_) * phi_u_plus(:,j);
        phi_f_minus(:,j) = (D_LL_ + exp(mu_sorted(j)) * D_LR_) * phi_u_minus(:,j);
    end

    % 
    % % ===== 能量流正定性校验 =====
    % % 校验图片中的公式：Re(iω Φ_fj* Φ_uj) > 0
    % fprintf('\n=== 能量流正定性校验 ===\n');
    % 
    % % 检查所有特征模式
    % for j = 1:2
    %     real_part_mu = real(mu_sorted(j));
    %     % 计算复共轭
    %     phi_f_plus_conj = conj(phi_f_plus(:,j));
    %     phi_f_minus_conj = conj(phi_f_plus(:,j));
    %     % 能量流计算
    %     energy_flow_plus = 1i * omega * (phi_f_plus_conj' * phi_u_plus(:,j));
    %     energy_flow_minus = 1i * omega * (phi_f_minus_conj' * phi_u_minus(:,j));
    %     real_part_energy_plus = real(energy_flow_plus);
    %     real_part_energy_minus = real(energy_flow_minus);
    % 
    %     if real_part_mu > 0
    %         fprintf(' ✓ (满足正定性)\n');
    %     elseif  real_part_energy_plus > 0 && real_part_energy_minus > 0
    %         fprintf(' ✓ (满足正定性)\n');
    %     else
    %         fprintf(' ✗ (违反正定性)\n');
    %     end
    % end

    % E
    E_matrix = diag([exp(-mu_sorted(1)), exp(-mu_sorted(2))]);
    % 计算E^N
    E_N = E_matrix^geometry.N;
    % 构建分块矩阵U（位移关系矩阵）
    U_block = [phi_u_plus,         phi_u_minus * E_N;
               phi_u_plus * E_N,   phi_u_minus];
    % 构建分块矩阵F（力关系矩阵）
    F_block = [phi_f_plus,         phi_f_minus * E_N;
               -phi_f_plus * E_N,  -phi_f_minus];
    % 计算柔度矩阵α_N = U * F^{-1}
    alpha_N = U_block / F_block;  % 等价于 U_block * inv(F_block)

    % ===== 提取子矩阵=====
    alpha_N11 = alpha_N(1:2, 1:2);  % 左上2×2子矩阵
    alpha_N21 = alpha_N(3:4, 1:2);  % 左下2×2子矩阵

    % ===== 计算透射率T_N=====
    % 提取矩阵元素
    A11 = alpha_N11(1,1);
    B11 = alpha_N21(1,1);

    % 计算透射率 T_N = |B11 + B12 * A21 * A11^{-1}|
    % T_N = abs(B11 + B12 * A21 / A11);
    T_N = alpha_N(3,2)/ alpha_N(1,1);

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
function D_beam = spectral_element_matrix(omega, material, geometry, L)
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




%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":45.5}
%---
%[output:1c1d99d2]
%   data: {"dataType":"warning","outputData":{"text":"警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  9.200181e-21。"}}
%---
%[output:84b53c01]
%   data: {"dataType":"warning","outputData":{"text":"警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  9.584400e-23。"}}
%---
%[output:4171ae8c]
%   data: {"dataType":"warning","outputData":{"text":"警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  2.162319e-18。"}}
%---
%[output:5651274e]
%   data: {"dataType":"warning","outputData":{"text":"警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  3.659351e-22。"}}
%---
%[output:2e8e965b]
%   data: {"dataType":"warning","outputData":{"text":"警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  9.193335e-17。"}}
%---
%[output:5d4f7ab4]
%   data: {"dataType":"warning","outputData":{"text":"警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  4.636191e-17。"}}
%---
%[output:878ee9df]
%   data: {"dataType":"warning","outputData":{"text":"警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  5.440218e-18。"}}
%---
%[output:7aa25792]
%   data: {"dataType":"warning","outputData":{"text":"警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  4.533684e-18。"}}
%---
%[output:8cfd2c89]
%   data: {"dataType":"warning","outputData":{"text":"警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  2.593166e-18。"}}
%---
%[output:5b85e6c2]
%   data: {"dataType":"warning","outputData":{"text":"警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  6.987622e-19。"}}
%---
%[output:2c863765]
%   data: {"dataType":"warning","outputData":{"text":"警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  5.139472e-18。"}}
%---
%[output:1dfd022c]
%   data: {"dataType":"warning","outputData":{"text":"警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  1.864015e-17。"}}
%---
%[output:75481f82]
%   data: {"dataType":"warning","outputData":{"text":"警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  3.525092e-17。"}}
%---
%[output:5664c51f]
%   data: {"dataType":"warning","outputData":{"text":"警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  4.975363e-17。"}}
%---
%[output:3c2354af]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAfoAAAExCAYAAACUOmA+AAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQmYFcW1PirKjLskhFWDC1HykrjFNyOiJioa\/IQYibIpiKhEYyQGZwAXFAUGGDWOOyJOIDqIiYnKc4kkGgkSJyb6cCM6UVFZMwlGYjJgRn3vr6Hu1DR97+2luruq7unvmw\/uvd1Vp\/+\/qv46p7YdPvvss8+IL0aAEWAEGAFGgBFwEoEdWOid5JVfihFgBBgBRoAREAiw0HNBYAQYAUaAEWAEHEaAhd5hcvnVGAFGgBFgBBgBFnouA4xAAQQ+\/fRTwl+nTp063LV161ZqamqiPn360O67775dCi0tLfSf\/\/yH9thjD9phhx18c8D0mI0bN9LmzZupb9++ee8LQtAHH3xAsKlbt26+6bS2ttJrr71GvXr1os9\/\/vNBktzunjj2AkPgBSy+9KUvBcofuHzyySe0zz77EOz\/5z\/\/SbvtthvtsssuRZ+X9wN\/cPfRRx+JZ\/y4konhmX\/\/+9\/bcfbGG2+IWw466CDx70477VQ0f76BETAJARZ6k9hgW1JBAKJ4yy230FlnnUX9+vUTIvC73\/2O\/vznP4v8ISivvPIKbdq0id555x0aNmwYTZkyhfbcc8+cfatXr6bLLruMZs6cKdLwXosXL6YVK1bQrFmzqLy8XHQWPvzwQ\/LOff3tb39L06ZNo9tvv52+\/OUvd0gGAgWh+sc\/\/kFTp06lNWvWCNE67bTThO1IV1733nsvvfvuu3T11Vdv1ynBPeh4TJ48mUaMGEGVlZW0YcMG+sUvfiE6B97r29\/+tujARLXXr2OD9164cCE9\/vjjdNttt1HXrl0Lcg3Rrauro5UrV4p\/\/\/a3v9HYsWPphhtuEPYXu9CpqK6upjlz5ohO1OzZs8UjkyZNyvsoMLnooovo+9\/\/Pp100km5+\/AsOhyw+de\/\/jVdddVV9NWvfrWYCfw7I2AMAiz0xlDBhqSFAMQNYiNF54ADDqCXX345J3rLli0Togpx\/cIXvkA77rhjzjR0EvD966+\/Lu6BlwcP+dprr6X9998\/d59X6NFpGDdunBCuoNfgwYNFRwEiDQE6+OCDhVf+8ccfCy8XHYT77ruPXnzxRWFLWVlZzls\/4ogjhNC99957ojMyffp0IXo777yz8OqPOuqo3Duj07Jq1SohbvgdQr\/33ntHthe2QJTXrVsX9FWppqZGdKjkBaFGGhBm4CA\/q0L\/r3\/9S3CjdnjU58MKPTh75JFHRKcL+MpLdhKqqqqE0OP9Ro8e7duhCvzCfCMjkCICLPQpgs1ZmYOAFPu\/\/\/3vNGjQIPrjH\/+YMw4iDu9YCh9+2HfffWnIkCHC+4dgw6uG8K5du5YWLFggvE54jvLy8+gRpke+MjKAz1dccQVdeOGFIpzd3NwsOg64ByLWvXt34dGjczFx4kSC0MArPu6448S9b731lhDlF154QQjUhAkTcqFpPPdf\/\/VfIiKhCv3XvvY1YfchhxwiOgXwtGtra4VYXnLJJbmwPyIQUe1FRAQYymjBX\/7yF7rrrrvo8ssvF17xvHnzhEd8zDHHCKFGBOCLX\/yiwBgXnrvuuuvo\/fffF7hCdL1CD4\/\/1ltvpZdeekl0YIAVQuzABRd4ufvuuwW26NiAD1xqZ0JigO\/x3A9+8AMh4Mcff7wQ8c6dOwu+b7zxRpEesHv11VcJZWbGjBmBhyDMKfVsSakiwEJfqszze4sxdIjM+vXrhbAjPPvoo4\/SL3\/5S4HOgAEDRMMPT1kKJ8aNIfQHHnigEBCI8HPPPUdnn302XX\/99b6owls988wzadGiRXTHHXcIwT7hhBNEFADe+emnn04PP\/ww9ezZky644AIhXAhVQ4AhxogGQOjxW319vYgowKY\/\/elPdOSRRwpRhJcJe2X0AR0ARCogkKrQy9A9DIWnOnfu3A42H3rooTR\/\/nzh0Ue1Vya4ZMkSuv\/++2nUqFEiL9jeu3dvMYSATgiiCBiugOcu50Cg4\/Hzn\/+crrzySoGzDLWrQl9RUSHuAQ549+9+97uCR+SB75HHli1bCJEKDEEAK+CDC78Bc3Q+0EHDMIAcVkBnD5EPdC5wDyIpsiyA\/xNPPFEM0+APoq8O5XB1YgRMRoCF3mR22LbUEMD47J133inE4bDDDhOeG7xEfI8QMDxuiIkMwXs9eoSUcb8qco2NjUKwunTpIjxmCAq+Q7gdwg9hwnfwzM844wwhavgNnQ48J8PHUuhPOeUUEfpHBwHepvRq4dkefvjhQughTvDuIVaf+9zn8obQx48fL0xFRwHChU7Om2++KcLWEHrYHNVepIu0MNyA99h1112FN68Kff\/+\/cXwAb6HB41xcYj98uXLxf8RFYCNfkKPIQHMmfje974nvHDZSVDH4QuN0UsOkTaEHhGR888\/n77xjW8ITDCMgwgE8MPEvyeeeELYhWEUTAYELsAewxx8MQI2IMBCbwNLbKNWBDA2j\/AxLogMRA4hcXjM8CKfeeYZ4S3\/6Ec\/oqefflo08BjzhveIRt7Po0cnAeKE0PR+++0nQsXqZDz5AngenQl4jIgIICwMccZz8FThScNTRBrykkJ\/zjnn0E9+8hM69dRTRfroPOA5iA\/yRlga\/6LTgdAyvGWE0CFK3lA2nkMUAe+PjgVswv8hllLokX8UexF6R\/7oqECI4Sn7CT3sRYcAgosJbgjf\/\/CHPxTvhTTg9XuFHp0SiC5EXnYOJE5RhR42AE+8K6IQ6NjJKAF4wpwNdFi+9a1viQmb\/\/u\/\/yuGBBBl4YsRsAEBFnobWGIbtSIADw5igb\/\/\/u\/\/pvPOO0809DfddFOHfHr06CG8UIwdQywRvs3n0d988830+9\/\/XszWRwgfQwCq0OP5P\/zhD2ISIIYKRo4cSc8\/\/7wYPjj22GNFBAF\/EHHMhkeIHR0KiJ4Ueozn\/+pXvxLL8TB0AFHC6oGlS5fSgw8+uJ1HD28Vz8I2jOvDW8ewAKIUWGb21FNP5RV6dDai2gsQ0XnChaVob7\/9Nj355JMiaoFIAbCBqMOjxwVspFcO0UUEAMMWGMqQQo\/vEeLHkAZC9ugkeJc8QugRskfnAnMTEPmQkyQR9cCF37Ca4NJLLxXiLWfwwwbwj3kJmOSIzhc6J1iNgQt5wX7MK0AHBGXDbxKg1oLKiTECmhBgodcEJCdjHwLeJVcQUExIg7cGzxrhXBm+heB85StfER71T3\/6UyHO6hg9Qr3oCEBIMCkOHrQq9BjLx9I3eM0I+2NZHSaLIeQOzxcT+XA\/BByijHsxzg5xgTDBe4TQY4Y9vG98j7F3RAQg9Oi0eD16iNizzz4rxskRRkdE4pprrhFhaSzPQ4dHevQYw4aISY8ewxRR7YW4oxODzkyxCxMekac6y10uBYTQYz4D3gG2I\/KAOQ4YwvBbwuc356BQ\/g0NDULoMfHwZz\/7meg0nXvuuaIjhjF4jO+DY0zIwzACohT4HrZgNQbKCq+pL8Yw\/24CAiz0JrDANmSCgFfoISjwuOH9wauEl6d6fxB6eM74DuIAcUUoF43+PffcI2bM43l4hQgxwyOW6+jzvSAEUYorvGjppSJ9iJmcFyCFHkMKEDyv0Pt59LAXY9knn3yymGWPkPTQoUPFGD7eEZPRIPSYTIhJavCyvaF7r91B7MUzCHnD+8YwAqIH6GhgXoG8kCfC9ZioiCiGekmhh8eMyYWYYIcOEYYhMOQivXDkgdny6BDIyXhIB52yKGP0eBZpo8PzzW9+M7cKAEKPYRxg9de\/\/lXMs0CnALzzxQjYgAALvQ0ssY2JICCFHiKKcVh4uwjpYnwWIXWEgTFRDuFfhL\/lOnmM5cpxW4T9EUKXG7NgWRzEGuPfqkcf1MNVX1TOgMd3sBFeP8QZtqBToXr08PIhcOpEMwxLYPkZltNBEBGqxtACQtJYTQCbMSnuscceE1EG2C2FHqHyIB65n70Iz8sLwghM4Rljtj06M8AGHQB0Qr7zne9s552rHj1WM2AYAFwAa3UdvYzAIA38RR2jx3wAzG2A3ehYoIOHPGXHRC5BRCQEnCCKghUNfDECtiDAQm8LU2ynVgTgMSMUC88P4eMHHnhAhN0RFkZjj1A2QvgIwSOce\/HFF+dmumPJFWa143mM1UJgpdBLI72T8eB9ym1Y1RdBKB5pwUPca6+9OrwjvFR8h53x5Dp6dDYgQPDKZfQB4g2hx0Q2hNuxbA3pqRv4QBQxdwDvjTFouRUs5gXA08fEPYTzMW6OfOPYq24whBfCsAUE+n\/+539Evpg8CFHG6ga\/ELwq9H6z7qVHj84LhlcQ3seKg6hCL21Epwgz7DF3AZ0c5A3RR6cO2GF7XExURKcF8xwQJcm3vbHWwsqJMQIxEWChjwkgP24fAhAerJdGA45wLELw2H0NIW6My6rjrhBIufUsJuVB3CGkEHaIojdELCeWoeOAyXlhQveqJ6yiqk7GgyePUDb+sEQPW99iTB9\/+A2dAnjR6AxgPBnvIoUKKwcQLsd4OLxULCOD6GPiH56Fx481734T3aQ9aug+n724F2KNIQx0QrCK4Te\/+Y0IscOjh6BiWR\/mIMAmfK\/uQe8n9Jhchw4X3hkdBFwI5aNDBU9bbnPr3RegUOmUY\/R4J0Q4wC\/KAezCsAHmUSA9dHrQoQCO6DxhgyVEVhBR4bX09tX\/UrSYhb4UWS\/xd0bDDm8Y48No3OEZYnkdxBm7nvldWNoGbxkihc4AGnh8huBilzcIAtbRQ3gQksaFiADC5N7Z4Wr6QYRTFXpM1kPoG4KDcWv8BhvgiUL4cEG04fViUiA6GuhwYF9+RCvQAZDb+CKagEgBdvyDjdiJDvYiXA7B87uC2Isd6rBiAB0OzBNAxASbzciDe9ABgNeM7WQx3o1hEgybQGzhIfsJPTpQEHSM0SPMjgvr2DEhDuP\/ctgBnTJ19zu\/d1CXIAIndObkToTSQ0ckB0M26PxhDgGiEBjOQH7YbAdRBCzxY4++xBsTS17fGqFHA4OdurwekjrTVvbQJfZqo+vdS9sSftjMBBCAaGBcFo13nAuChD94o\/JENRnyxudCJ6XJfHE\/RAYhem\/IW96DcDvC\/pggCEFTJ+oVsx\/vCa\/e7\/Q9v++BDezIZ0sQexEtgJhKewvZ6H033Cufhw1BMJTpBzmhDvfKA4aCnoSn2h\/k\/Ytxwr8zAmkjYIXQy4lM8pAPuX5V9S7Q21ZnDKshVYCqHnCRNsicHyPACDACjAAjkBUCxgs9xBuzmrE9JbwE1aNXJ994j+H0TobCvQi5FgvrZUUE58sIMAKMACPACCSBgPFCD68d42he4ZbCjuVBEG\/vZ+8a6SDnUScBMKfJCDACjAAjwAhkiYDxQi\/BySf03tO4pNfu9eDxvFxr7Ac4JvHIE66yJITzZgQYAUaAETAHAWwmhT+bLxZ6IiHwmEWLGcB8MQKMACPACDACEgEcNoWzF2wWe+uFXkfoXk72A5nYDIOv9BFAJwsbxzAH6WOv5sg8ZIu\/zJ15MIsH74ouM6wLboW1Qo9XVMPzfpPx1FB9ocl4UuhtJzM47ebdyRyYwQnzwDyYgYAZVrhSH6wWel3L61wh04yqEc0KdMqw29iYMWPE6gi+skGAecgGd2+uzIMZPLiiDVYLvfTq5baXUTfMcYVMM6pGNCtwaAmOB8Xe8WVlZdES4adiI8A8xIZQSwLMgxYYYyfiijZYI\/SxGSuQgCtkJolR0mlzw5Y0wsHSZx6C4ZT0XcxD0ggHS98VbWChJyJXyAxWdM28ixs2M3hhHpgHMxAwwwpXtIGFnoXeiBrFAmMEDeLcdx5CyZ4L5iF7DmABC70ZPGixwhUytYCRUSLcsGUEvCdb5oF5MAMBM6xwRRvYo3eo12ZG1YhmBQtMNNx0P8U86EY0WnrMQzTcdD\/FQq8b0QzTc4XMDCGMnTU3bLEh1JIA86AFxtiJMA+xIdSSgCvawB49e\/RaKkTcRLhhi4ugnueZBz04xk2FeYiLoJ7nWej14GhEKq6QaQSYEY3ghi0icJofYx40AxoxOeYhInCaH3NFG9ijZ49ec9WIlhw3bNFw0\/0U86Ab0WjpMQ\/RcNP9FAu9bkQzTM8VMjOEMHbW3LDFhlBLAsyDFhhjJ8I8xIZQSwKuaAN79OzRa6kQcRPhhi0ugnqeZx704Bg3FeYhLoJ6nmeh14OjEam4QqYRYEY0ghu2iMBpfox50AxoxOSYh4jAaX7MFW1gj549es1VI1py3LBFw033U8yDbkSjpcc8RMNN91Ms9LoRzTA9V8jMEMLYWXPDFhtCLQkwD1pgjJ0I8xAbQi0JuKIN7NGzR6+lQsRNhBu2uAjqeZ550INj3FSYh7gI6nmehV4Pjkak4gqZRoAZ0Qhu2CICp\/kx5kEzoBGTYx4iAqf5MVe0gT169ug1V41oyXHDFg033U8xD7oRjZYe8xANN91PsdDrRjRierNnz6a5c+eKp8ePH0+TJk3KpbR48WKaMmWK+FxTU0PDhg3zzcUVMiNCaMRj3LAZQQMfU2sGDcyDITy4og1We\/QQ8hUrVtCsWbOopaWFxo0bR8OHDxeC3tTURNXV1TRnzhxRZOT\/+\/btu10RcoVMQ+pGJDNY6CPBpv0h5kE7pJESZB4iwab9IVe0wWqhhzePS3rx6me1E1BeXk74rU+fPr5evStkai\/lKSbIDVuKYBfIinlgHsxAwAwrXNEGq4Xez6OH6FdWVgphz9cJ8BYhV8g0o2pEs4IFJhpuup9iHnQjGi095iEabrqfckUbrBZ6kCqJ6NmzJ9XX15MMzXs9eHQKVq9e3WEMXxYKmcaECROooqKCunfvLv74Sg8BNGwbNmygbt26ESIwfGWDgGs8lJeXCSBbWrZkA2jEXF3jISIMmT2Gtgh\/a9eupaqqKmpoaBAOpK2X1UIPMV+3bp0Yo8c1efJk6t+\/vwjPRxF6SeKYMWNo9OjRtnJqpd1bt26l5uZm6tq1K3Xu3NnKd3DBaNd4OOCA\/QUtQ4d+RLW1zdZQ5BoP1gC\/zdCFCxfSggULcmaz0GfE4KZNm8TkOxmql949BH7+\/Pk0b948YZnf+L3XZOnR19bWUq9evdijz4BTTKbcuHGjwL6srM0L4yt9BFzjQXr0o0a10j33tKYPaMQcXeMhIgyZPSY9+sbGRqqrq2OPPismign90qVLO4TqeTJeVkwFy5fHJIPhlPRdLvGwejXR\/m0OPZ17LlF9fdLo6UvfJR70oZJ+SjxGnz7m2+XoF7rHWD28eF5eZwBBIUzghi0EWAne6hIPqtB\/4xtEzzyTIHCak3aJB83QpJocC32qcPtnhvAWxuWXLFkibuANcwwgJaIJ3LBFBE7zYy7xwEKvuXCUYHIs9A6RbhKZzz9PdPTRbeB+9plDIBd5FZcExmbWXOKBhd7mkmiG7SZpQxxErJ51H+fF1WdNIvOb3yT67W9Z6HVxy+mEQ4CFPhxeSd3tEg9JYZRGuiZpQ5z3ZaE37FAbFvoeoWbd77BDaUU+4lT2IM+6JDDs0QdhnO8phAALvUPlwyQyWeiDC73Eqk8fonfecahAZvgqLgk9ImMoI7h4Ml6GhcrirE3ShjgwskfPHn2c8qPt2SgCw0KvDf5cQlF40G+FnhRZ6PXgWMqpsNA7xL5JZLJHzx59llXLVaG3LerjEg9Zlue4eZukDXHehT16wzx6bPCBsUVcPOu+cNEuVaziVPhiz7okMKpHz0JfjHn+3Q8BFnqHyoVJZKrihXFnNFClcEURGBZ6\/SUjCg\/6rdCTIgu9HhxLORWTtCEOD+zRG+zRs9CzRx+nckd5loU+Cmr6n3GJB\/3opJciC316WCeek0lkskcffIyePXr9VcMlgWGPXn\/5KLUUTdKGONizR88efZzyo+3ZKALDQq8N\/lxCUXjQb4WeFH\/yE6KxY9vS4jF6PZiWWios9A4xbhKZ2ABGXhy659B92tWMhT5txP3zc4kHMxCNZoVJ2hDtDdqeYo\/eMI+ehT546F7FqpRWKMSp8MWedUlgVI8e721TGXGJh2JlzuTfWehNZiekbSaRyULPQh+y+Gq93SWBYaHXWjRKMjGTtCEOAezRG+zR4\/xsbN1ZClcUgSnVTlGS5SEKD0naEyftadOIrr22PQX26OOgWZrPstA7xLtJZKrixUJfuJCx0OuvhMuXb6Vjj+1MLS1bQh0upN+S+Cmy0MfHsNRTMEkb4nDBHr1BHr162hZIZaFnoY9TuaM8KztPixa10vDhnaIkYcwzLPTGUGGtISz0hlC3ePFimjJlirBm8ODBNGvWLCovLxef1d9qampo2LBhvlabQiYL\/Xrq0SPYGL0Xq1JaoZBU1VMxPfdcovr6pHJKJ10srcM4vbw4dJ8O7i7lYoo2xMXUao8eJFx++eVUX19Pffv2pdmzZws8Jk2aRE1NTVRdXU1z5swR38n\/4z7vZQqZLPQs9HErdJznXRd6mzqDLs2ViFMms37WFG2Ii4PVQg9h79Onj6+nDm9+xYoVOQ+\/0L2mkMlCz0Ift0LHeZ6FPg56ep9lodeLZ9TUTNGGqPbL56wV+k2bNtG4ceOE915ZWbkdDqp3jx+9n9UHTCFT3bIT9vEYff7izaH7uFV\/++ddE3r1yGe8LXv0+suM6ymaog1xcbZa6CdOnEhnnnkmzZw5k9atW9dhjN7rwcPDX716tegY5AvdT5gwgSoqKqh79+7iL+1r+fJONHBg+wSoefNa6eyzW9M2I5P84MFs2LCBunXrlptjUciQNWs6Ud++7VitWrWlZE76S4ogFdNRo1rpnnvsLnuDBpUROs\/yeuKJLdYsVw1bH5IqE6WaLtoi\/K1du5aqqqqooaHB16G0BR+rhR4e\/X777SfC87gmT55MPXv2FGIeReglaWPGjKHRo0enzuHzz5fRyJE9cvnW1jbT0KEfpW5HFhlu3bqVmpubqWvXrtS5c+eiJkCUjjtu39x9DQ3rqbJyS9Hn+Ib8CKiYotyh\/Nl8jRjRgxoby6ysT2Hrg808mWj7woULacGCBUr7wkKfCU9+oXuEWSDw8+fPp3nz5gm7pAcfJHRfW1tLvXr1ysyjh\/cBL0RepeTRt7S00MaNGwX2ZWXtGOQrXF6sli5tpQED7PZAM6lISqYqpi549P36lRGGI2ysT2HrQ9Zlx7X8pUff2NhIdXV17NFnRTAqAjz4ESNG5EIqEPpFixYJD\/\/RRx\/tEKq3YTKed8tO7Op1zTVZIZxuvmEnH5XyfIakmFExdWF5nXq6ITCzac5L2PqQVJko9XR5jN6AEqDOrIc5EP7+\/fuLWfg2Lq9joQ8+656FXn8FVDHF1ssQRpsvdedEvIdNHWcWejNKHgu9GTx02BRn\/PjxHSbb2bZhjncnLxe8qqDFJGzD5hV6m2ZUB8Uk7ftcEnrvqgxgaVN9Clsf0i4rpZIfC71DTJtCJgt9cI\/eG\/3ALm5oyPmKjgALfXTsdD\/JQq8b0WjpmaIN0axvf8raWfdxX1x93hQyWehZ6HWW67BpuST03ogPsLBpOIKFPmzpTeZ+U7Qh7tux0Bt0qI13b26bGqa4BTFsw+btFLFHH5cBEmvOscmMbaLo9+Z+Qt+nT9umOTZcYeuDDe9ko40s9DaylsdmU8j07uTFQp+\/kHmF3qaJVqZWHZeEXi0fEHi5zM6Wg21Y6M2oJaZoQ1w02KM3yKP3Cr1NHkjcghi2Ybv3XqJx49pzTVvo5YxuW4QjCD\/qvAfbO5ned5E75NkyaTNsfQjCL98THgEW+vCYGfuEKWR61\/0CMJeEpFABCNuweYc5kpxRLXn5\/veJbrut7ehT5I\/LprXZxSqgS0KvdpoxrCP5YqEvVgr4dxUBU7QhLivs0Rvk0Usv8QtfIPrrX9uotaVhilsQwwq9N\/px3nlE8+fHtcL\/ecmL9HLVsLCrQm97NEmWD7wHNp2SQm\/LXI6w9SGZks+pstA7VAZMINNv3S8Lff5C5hX6JD16Fnr7KrvKGcQdURlcv\/410Yknmv8+LPRmcGSCNuhAgj16Az16dfKQLR5I3MIYtmFLc+IiC31cdtN9Xu00y7kbksMkO4Q63zJsfdCZN6fVjgALvUOlwQQyvTOe5eQhFnr\/guadz5BUqFkVDZkHh+7NrvxqXZJDX2p5sWHeCwu9GWXMBG3QgQR79IZ49Kp4qJOHbPFA4hbGsA2bdx9z5J9EA85CH5fZ9J\/364ipkzdtmPcStj6kj3Jp5MhC7xDPJpB5991E48e3gYqGSI4pstBvX9DyzWfIQuhdirh4157bsrmMt4SoE\/HkO6hevg2csdCbITAmaIMOJGJ79PJc+JUrV4ay59BDDxXnxnfp0iXUc0ncbAKZ3sZJhhqTCkkngWOcNMM0bGlOXFQFwi90b4NoBOXFBaFXy4a6F4DK4\/XXE111VVBUsrkvTH3IxsLSyNUEbdCBtDahnzRpUu5c+GKGATycD89C346UFHbZOKmhxiQ81WIcpf17mIZNbbQx2Qp\/uJIQXRb6tEtCvPy8ZQNL6+Rl0zh9mPoQDzF+uhACLPTb0JEePQt99AqjeiEjRhA1NHTclMWGMcXob9\/2ZJiGTfU8VaFPYk27N3qATpd3PoUrp+Z5txW2sYNZaKKkTeP0YepD3LrHz+dHgIXeodKRNZl+y4HUXcqS8FRNoy9Mw4bNcYAJrqTnM3gPR2GhN63kdLTHb3xe3mHTOH2Y+mA2I3Zbl7U26EIvduheGoJQ\/Ny5c8XHmpoaOuKII2js2LG0bt068d348eMJXr+JV9Zk5vNCbFv7G4fbMA2btzFPEifvufcuC713W2HbPPp84\/NquZRlxfS5L2HqQ5x6x88WRiBrbdDFjxahh8hD0GfNmiXsmjx5Mi1ZsoQaGhrEuH1LS4v4rn\/\/\/jRs2DBdtmtLJ2sy83khNo0pxiUjTMPm3cAmyYmLfkKvCqJL0Rav0Ns2ZBTkDAKM2V93XVtpTWKoJ249kM+HqQ+68uR0tkcga23QxUlsoZciPmLEiNxkPL\/Jdvhu0aJFojNQXl5Lxd0kAAAgAElEQVSuy\/5cOn55Ll68mKZMmZKLMuTrZGRJZiEvBKZv6zs5v+d90IbNb5gjyYmLfuPWLPTaq6+WBAuF7WUGavje5KWrQeuDFuA4kbwIZKkNOmmJLfR+k\/HyCX1SM+2lDQBGzuRvamqi6upqmjNnjsBL\/r9v377b4ZclmX7C5dcoueQ5+hXgoA2bipf0yJ5+un3\/ct1eaCkJvXdbYd1Y6my4vGmpAv6jHxHdeGP+3NT3NPUdg9aHJDHltImy1Aad+Dsh9PDcH3\/8cdq8eXNO6PHdihUrchEEdDL69OnjO3SQJZknnNAWQsTl1+h4w9Q6yTcpraANm9pIS6FPcuJiKQu9yaFtb9n1Kxf5yrfaKRg3juiee0yqCW22BK0P5lnulkVZaoNOJK0XenjuCxYsoBNOOIFuvfXWnNBD2HHJCYDezyqIkswJEyZQRUUFde\/eXfylcZWXl4lsevdupaam1u2yvOSSnWj+\/J3F96tWbSFMInLxQsO2YcMG6tatW8GhnUGDyggNtRcvieOoUa10zz3b4xgVs+nTO9GMGZ1yj4MDfHf\/\/W3fzZvXSmefrS+\/qHbqeE5iK9N64okthH0dTL\/WrOlEffu28TFgQCstXVqcj4EDO9Hy5W3P4H48Z9IVtD6YZLNLtqAtwt\/atWupqqoqN9\/M1ne0Xugh4Mcff7zAXx0a8Hrw8PBXr17tO\/NfCr0kccyYMTR69OjEOX3++TIaObKHyKe2tpmGDv1ouzwfemh3qqrqKr7H77jPxWvr1q3U3NxMXbt2pc6dO\/u+Ihr0447b1xcLfI\/fcb399jvaIAL24EBey5a9T3V1++S+y8ebNgNSTEjFENk2NKynysotKVoQLasRI3pQY2NbhzmozWrdQ6cRvJp0BakPJtnrmi0LFy4UDqS85MRyW99Tm9AH2QJX97a3EOhnn31WiLd3XkAUoa+traVevXql5tFLrwJeOjzFfFe\/fmWEsWlcLS3mN7xRKgMmdW7cuFFgX1bW1mh7L3jy8Dpxeb2w889v97J1Rj7UdJGvyx69Ws7wrjZEK9QygegDohBBr4svLsvtx3Ddda1UVWWOVx+kPgR9T74vPALSo29sbPz\/jn0de\/ThIdTzBCrCjBkzCN43Jtj5CT1yChO6T7PX5j2WVo7T5xM4jEHismncNAzTQcYks9gMxW\/JGcbtMS8AlzzvPMy7mnqv90RA09\/Ne7IgJqyGHWpQ39mkiXlB6oOp5cglu3iMPmM2MTavbsgjzenZsyfV19fTiy++2CFUb9pkvLAzf9UGybaNTIIUlWINm9oxuuACIpz2572S2DjHOxMdHS1E9EpB6E1efgbuv\/Utol\/9Kl4H2LvzoSliX6w+BKlTfE98BFjot2EY5vQ63aF7lUavR2\/y8jp1ljiEBMvDil1qgwShSWEKQTGTtP5erGEL0jEKck9Yo\/2EHh49+HDJo\/c7EdBkodd50p53UyQTombF6kPYcsz3R0OAhT4Pbpj09sADD3Q4mS7KwTdhabFlwxw0qBAP\/Iux+TDhRnWnPFM8j7A85bu\/WMMWZOvSIDujhbVXxRzPgi90tEpB6E3dJlblOWwdysf\/smVE2+b0iluyrl\/F6kPYcsz3R0OAhd4Ht0KCbtrRtN5owMiRI1OZcIGzsKdObcsdO9\/NnBmuAAYRvHApmnF3oYYtzBpp3fsOeMetXRV6NWKEWehyBYNpw0TeyINO77uujuiHP2yvDzrTDlvLWOjDIpbM\/Sz0LPShS5bXE4HXEPZSRS9o2D9sHlncn69h8064KobZ1VcTTZ\/e9gY6Gmq\/CWoQRdc8elXo99rrU\/rwwx2N8GzVsugV+SR2i3zjDaJDDmnP9fTTiX75y\/RrBAt9+pj75chCn4eHQqH74cOHl+yhNjpmCEvI1XCyDjEzoUr5NWzqMAdsDBJO1bmXud+4tasevdoJxX4Ncu+AJMQ0SnnzcvHjH3f0vqOkme8Z7wQ9DA+gnqW5WRULvU5Go6fFQl8AO+8GNLg1zaVrYWlNmkyvYMUV5yTDl2Gx03W\/t2HDO158MdZFt+UQZpgjTKi\/kP3eBh\/3YoIabHPNo1cnt2ETILlJkwlC750sh4Oekj7x2ltnwf0ZZ7TtoZ+G4LPQ62pZ4qWTtDbEsy7407E3zAmelbl3Jkmmt8G46ioijNPHvbwiBEHEciNbL2\/D9uST2Byn7W3CTgoLs0dBXKHH86aNY0cpA9XV2J2x7UnsEid3IMxy5j3qDpbPfe977W\/01FNEAwdGecNoz9xyC9GECR2f\/elPic4+O1p6QZ9ioQ+KVLL3JakNyVreMXUW+v9fj5sUmV4x1n2AhtezP+00oiVL0iw++vJSG7aFC8to\/Ph2kY8SNtXh1XsPtIFFckMW6dHjOxO83rhMqHhhC2F1O9wsOjLeDjI6ezhLHh2PtC\/Ygo2TVM5hwyOPEH3ta8l4+Cz0abPsn19S2pD228UW+ihL50ybgZ8EmerRqSB1xgyiK67QT69X7LMYT9TxVmjYXnihmZ58sivNnNm2zW2cd1E7WWEjAvJ9VPFDGnJJJP7vqtDLfd9\/9asedNFFbTzEHWoKUz6AsRR5+ZyuJXRh7PC7F0MI6PzJ7ahV+6J0RgvZw0Ifly09zyehDXosC5cKC71mj95vbC\/phtIvz4ULic45J1xhyPJuNGw4gUwu64oqzuo7qB45PPFC2wz7vbu65a4q7kjLNaGXEzwrKrbQokXr6Y03euTOFUDo+uab0ykd3g2KdJQD3ZaDe3j4foKPXRtHjozv5bPQ62YtWnos9NtwC7Mzngp1krvkhaVUB5mo9IsWdfTa0\/ZE7r2XCMMDSXoaYbEtdn8+D06XhxQnhK+uyQeXcttb7zuZvid8EA4g9LgGDtxCc+eupx49epB6yE2QFQ\/F8sn3O8qAen4A7ku77kSx3c9ute4NG9Y2vyDK5D0W+iiM6H9Ghzbotyp8irE9+vBZmvdEHDJR2f\/0J6Lvfrfje8Hrw9htlEoeByE\/7x42zJ9PdMAB6dtT6F38QqG6cYu6l7n6HDodzz7bdoCN3+WS0E+atIXGj28TehzlKg9TmjiR6IYb4pRM\/2e9hwbhLozDYzw+7boT9e1kZ9XPy1eFH5Mdv\/71YO\/FQh+VDb3PxdEGvZbES42FPkLoXobs\/CbomOKJ+C0NM8VT+vnPic48s2PBxdjwnDnN9J3v7JP3mNqoRT3KckRvJAB5S9Hz2mG70Ktl5c47t9App7QJPY4LTmLPBuSnnhegiqGuSE7UsqLjuXyhfTVttBMnnkiEVTiyXqq\/s9DrYCJ+Giz08TE0JoWgZOZroEwRUC+gfmFxb2Mjx62T9p7uuIPoZz\/bfuYy8sUSqv3220Lr17cLjO7CEUbs1XvlGHG+jhPstF3o1bkMOM\/94IPbeQiDWyHOkM4rrxANGbL9XXEmXeouJzrTk\/VPPe2wUPrAAVG3qirctYV2262Zjjqqq\/aOr853NDUt6YzFbdeCaoOpOEi72KPf5tGfdVY11dbOoW7dKkXIcOVKoocf7rg5ipdM3CeX\/YQ9BzvtglGok6J6VFjWVlnZtoxMzjIPYyvyASa33kr04ovbC7ual3qgTxoejHfjFcxpQFTGe115ZfsZBOrSOe92uPK5LNeah+Em371q+LylZfsOl1fsv\/MdoptuKh6CxnPYMwI4+10oYwjRm153dGCMNKTwL19OhK2ag15SrPAvnkNnAGmVCm4qTuoESMwbwXLil15qu8O7\/BHfxV36ykIftJRacB\/IPProykCW2iTufi8kGxu\/0GkhAPL1jL0zj4ulgZA9drzzppeG0MvGVk48w2fpTcr\/o2E477y2t\/DO+PaeZOeK0KurC1at8o+seDtJ8t1lfZCfZfnKVw5wP4QfWMb1tgJVWMNvkvVHLtvzE6sgr6B2BnD\/yScTHXYYUb9+bU+rHYMoHfggNgS5R\/W0ZVmB7Z9+2namwKpVRG+91W5zmPbFL38W+jZU2KPf5tH7Cb2sPOg5jxnjXg9aViI5ezhqI1OoUQ+KXVpCL23NJ9qqgHkP0PGbOIb7bffo1dUFCN0XGkLxLn8L0rijHmGHua9+lcU9CF6oj6gPzc3N1NjYnR57bOdcNCDI83HvidMBiyvMUWxXOzn4\/9ChRF\/5SntKcd6HPfoojBj6DMg844w6mjOnlgYM6J3z5gw1N1GzVPFHBfnNb4hWrCB6++382eI+nPhVUdHekIetXGkLPd6m0IRFv1Py8nm1Ngu9GpaH9zN8ePG5EmpUSPXKZAlB5w5HMb\/7rnud40Qrn5J4odMcUbdkp\/xf\/yL6xS\/a66fuznpa76vmo7YdMmKE4Qq\/fUHCtjNh34eFPixiBt\/vCpkGQ1zUtCyEXhoFsZINJBqOQmOffuP0Ngu92tkJKvRFyeQbYiOgoz7IEL3aGZMd+S9+sW1yZFMT0fvvE334IRHKdrFjoIO+GOrRTjsRde9OdOSRRHvsQdSpU5sj4DdRTn6ftHAHtV\/e54o2WB26b2pqorFjx9K6desEL+PHj6dJyrFWODJ3Co49I6Kampq8R+S6QmbYQmzS\/ToatrTeB94FGkrZObBZ6NUoBVZgVFYW9+jTwrmU87GpPrjMkyvakIjQF9stT8eueN499uVneeY9OgHV1dU0Z84cUQ7l\/\/v27btduXSFTJsrnG0Nm3esOouDX3TwrU7EgzdnGw86MDAxDebBDFZc0YZEhL4QRbNnzxY\/q563LkrVtOHNr1ixgmbNmkXl5eWE3\/r06ePr1btCpi4cs0jHtobNOzEvyS1ik+RDTkqUKwxs4yFJbLJMm3nIEv32vF3RhtSEXobZb7jhBqrEQu0ELlXovR2KQh0MSeaECROooqKCunfvLv74Sg8BNGwbNmygbt26iY6ZDdddd+1Ml122kzAVy9JMG18MgmF5edsJdQMGtNLSpa3Co7eNhyDvads9zEO2jKEO4G\/t2rVUVVVFDQ0NielWGm+aitBDZDGOLr3rJF7M25HwevDw8FevXu0bSZBCL+0aM2YMjR49OgkzOc08CGzdulUsJ+ratSt17tzZCpweemh3qqrqKmxdtux9wja+Nl3Yy37kyB7C5IaG9WJ83kYebMI8qK3MQ1Ckkrlv4cKFtABbGm67WOgL4JyGF4\/s5fg8IgVySCCK0NfW1lKvXr3Yo0+m7hRMtaWlhTZu3Ciwxx7rNlzTp3eiGTM6CVOxo5xt1yOPdKLhw9vshzcPr95GHmzDPYi9zEMQlJK7R3r0jY2NVFdXxx59PqjT8OLziTy+jxK6t73XllyxTz5lG8ck1T3ibZyM552IB5Zt5CH50pl+DsxD+pj75chj9AV4SGPWvSrycqa9apI3VM+T8cyoOPmssLFhs1no1Y1ysG+APNzIRh7MLtnRrGMeouGm+ykW+m2Iepe56QY6X3oIbU2ePJl69uzpO+7Oy+vSYkJPPjY2bK4IvXr6no086CmBZqXCPJjBBwt9xkLv3SxHFovBgwfnJv3xhjlmVJYgVtjYsNks9Krt6tJAG3kIUr5su4d5MIMxFvqMhV5nMXCFTJ2YpJ2WjQ2bzULvNz4Pzm3kIe2ymkZ+zEMaKBfPwxVtiL28LqvQfXGKgt\/hCpnB39i8O21s2GwWevXEOjk+z0JvTr2wsT6Yg54+S1zRBm1Cv3LlyqLo6tj6tmgmEW5whcwIr27MIzY2bLYKvToRb\/FiorPOai8GNvJgTCHWaAjzoBHMGEm5og3ahB7r15Pa8S4GT4EedYXMQC9r6E02Nmy2Cr16Yh28efW0Pht5MLRIxzKLeYgFn7aHXdEGFnoicoVMbaU7g4RsbNjyCT3Ozb7vPiJT19arB\/J4bbSRhwyKa+JZMg+JQxwoA1e0gYWehT5QgU\/6JhsbNj+hz7c+PWn8wqTvPchGfdZGHsK8uy33Mg9mMMVCv40HnoxnRoG03QobGzY\/oVfPd1c3ojGJHzkR79xzierrO1pmIw8mYavLFuZBF5Lx0mGhj4efUU+7QqZRoIY0xsaGzU\/o1bC4iUJfaHwelNnIQ8iiZsXtzIMZNLmiDbFD92bQEc8KV8iMh0K2T9vYsEmhl2e5q2F7oGmi0A8aRPTkk21ceyfisdBnWwd4CMUc\/KUlrmgDCz2P0RtRu1jo06Eh30Y5MncbeUgHuXRzYR7SxTtfbiz0ZvCgxQpXyNQCRkaJ2NiweT16NZRvokcfZKKgjTxkVGQTzZZ5SBTewIm7og3s0bNHH7jQJ3mjjQ2bV+jV8XnThV49yIZDxkmW7Ghp21gfor2p2U+x0JvNTyjrXCEz1EsbdrONDZsq9BjvxrI19TJtjF6NOPiNz8N2G3kwrChrMYd50AJj7ERc0Qb26Nmjj10ZdCRgY8OmCj2WqcGjx4XJeQiTmyb0xcbnWeh1lGQ9adhYH\/S8uVmpsNCbxUcsa1whMxYIGT9sY8OmCj3WpCMcbrJHn+8gG9VmG3nIuOgmkj3zkAisoRN1RRvYo2ePPnThT+IBGxs2VejhvWOzHFM9enUi3s9\/TjR0qD+LNvKQRHnMOk3mIWsG2vJnod\/GA4BYtGgRzZgxg3bbbTfaYZvb8M9\/\/pP+\/Oc\/0\/7770\/XX389XXPNNdSlSxcz2PNY4QqZRoIb0CgbGzZV6BGux2Y0+BeXaaF7Vejzjc\/Dbht5CFjErLqNeTCDLle0IbZHL4X+6KOPptbWVho+fLhg6Pbbb6dPP\/2URo0aJToBWQj94sWLacqUKcKempoaGjZsmG\/pcYVMM6pGNCtsbNhUoYeQ4oJnj\/+bJvTqRLx33mnvkHjZspGHaCXO7KeYBzP4cUUbtAl9dXU13XTTTTRkyBDac8896d5776WZM2fSf\/7zH5o2bVrqQt\/U1ESwac6cOaLEyP\/37dt3uxLkCplmVI1oVtjYsHnXzePNMVYPzx5CL3fMi4aI3qeCTMRjj14v5nFSs7E+xHlfU591RRu0Cf2sWbPoww8\/pFWrVtFHH31Ee++9twjjH3LIIZmE7uHNr1ixgmBXeXk5zZ49m\/r06ePr1btCpqmVJYhdNjZsfkKPCXkYqzdN6INMxGOhD1JS07nHxvqQDjLp5uKKNmgTeoTue\/XqRccccwz98Y9\/FB78hRdeSPgeYfO0Q\/cQdlyTJk0S\/3o\/q8XFFTLTrQJ6c7OxYfMTeoTF4T2bJPTq+LzfiXUqkzbyoLckmpEa82AGD65oQ2yh\/93vfkcPPvggXXzxxSJMfsYZZwjP\/tVXX6Wzzz6bevbsmUno3uvBw8NfvXp1Tvj9hH7ChAlUUVFB3bt3F398pYcAGrYNGzZQt27dRATGhmv69E40Y0anDqY+8cQWuuiiMiH0vXu3UlNTq7ZXqa\/fmR58cCeRtrzOOusTGjjwMxowIH8+y5d3ooED2+xcurS14L028qANYIMSYh6yJQNtEf7Wrl1LVVVV1NDQQJWVldkaFSP3WEKPl8f4O8QRIv\/ss8+KcXlMxHvzzTcJM++feuopWrNmDfXu3Zt22WUXMX5\/7bXXitn4SV5RhF7aM2bMGBo9enSS5nHaHgS2bt1Kzc3N1LVrV+rcubMV+NTV7U11dft0sHXZsvdp5MgetGZNJyH0+Bz1kmmMGNGDGhvLCiZTKC\/VzoaG9VRZuSVvWjbyEBVfk59jHrJlZ+HChbRgwYKcESUt9BiP\/8tf\/kJXXHEF7bffftS\/f3\/aa6+9aN26dYSCetRRR4lOQG1tregV4TeM2++xxx7UqVNHT0g3rVFC97ATww\/s0etmo3h6LS0ttHHjRoF9WVlhUSueWjp3+Hn0LS1bqF+\/No8ek\/FWrcovqoWsxPPwwiH26oU05RI+TPpTL4g9PHb5u\/ztxht3oquu2ll8hD3e39U0bOQhHbbTzYV5SBdvb27So29sbPz\/znxdaXv0AEcur\/v+978vPPUjjzySzjnnHLr00kvp9NNPp5NPPjmT0L03VM+T8bKtOMVyt3FM0jtGL2fZI1gVZ4weAi630wVuSHfkSAh\/2\/I99UI+F16IkHzbt34z\/YPOuMfzNvJQrGzZ+DvzYAZrPEa\/jQcp9Jjdjl4oZtzDM0PYHtdhhx2WidDz8jozKkpQK2xs2PIJfRhh9Yo2ooXqVrpjxxJNnZp\/3bt8fsgQoiVL2j7NmoVJqO0ph+l42MhD0DJm033MgxlssdB7hB5j848\/\/jh98MEHdNBBB9HKlSvpBz\/4AW3evDkToYd5vGGOGZUliBU2Nmz5zp+PIvTwzCHqMhwPzxwH5Xg9+EJYSkHHPeqmOEGX1rFHH6SkpnOPjfUhHWTSzYWF3iP05557Lt11111izfyuu+4qwvjjxo0Ts6iz2DAnTHFwhcww72zavTY2bF6hf+QRInjWYYUeIn\/99UT33tvGCkQe29QWGkv3408N+ctldOrSunxn0Ktp2ciDaWVZhz3Mgw4U46fhijbEmnUPGGXoftCgQdSjRw869NBDBbpLliyhV155hS644AKaPn166uvow1DsCplh3tm0e21s2LxCL\/eQDyv0d91FdNFF8URe8ql69Z991hYhkOP9LPSmlfr89thYH+xBN7ilrmhDbKH\/+OOPxQQezKSXB9oARozX77jjjgJRHG6DHfJMXTblCpnBi695d9rYsOkQetXjjurJq2wi\/C9P0UP4HulLoS+0x71Mw0YezCvN8S1iHuJjqCMFV7QhttBLMHGgzb\/\/\/W\/afffdcwKvA+g00nCFzDSwSioPGxs2CCqEVV5RPHo5fh5lTN6PC9UmjPG\/+2775L5Cp9ax0CdVsqOla2N9iPamZj\/lijZoE3rMcseEvBtvvNHY42jzFSlXyDS7yhS2zsaGzSv00mMOGrqX9wGZICIclF\/ZecA4PToQchY\/e\/RBEcz+PhvrQ\/ao6bfAFW2IJfQvvPACLV++XKCL2fbYGe+UU07Ju4UpNsnBNoJf\/\/rXO4T59dMTLkVXyAz31mbdbWPDFkfovSF7iLCuSx2nR2cCnQhcGLMvdtnIQ7F3svF35sEM1lzRhlSFHkfWYg98nE+\/7777msHktgmFI0eOtH73I2MAjWCIjQ1bHKFPypsH9JddRnTzzW0kwKMPs3mPjTxEKG7GP8I8mEERC72HB7\/QPSbkvfXWW4Qz4DERD58h8thL3u9c+KyodYXMrPDTka+NDVtUoVdnwmOdvPS4deCINNRJgiz0ulBNNx0b60O6CKWTmyvaEMujV6FWhR4H1zz99NPioBvsf3\/ZZZfRPvvsI4T+vvvuo8GDBxt1OpwrZKZT9JPJxcaGLarQJ+nNgx3vFrr4LmiHwkYekimR2abKPGSLv8zdFW1IROgffvhhevfdd+miiy4yStDzFR1XyDSjakSzwsaGLYrQq2PzQcU3LKIs9GERM+9+G+uDeSjGt8gVbdAm9O+88w7dfffdVF1dnTulTl1XHx\/y5FJwhczkEEo+ZRsbNq\/Qy8luhWbdY2McbJCDS+dMey9Dcua9\/F7ulFeMSRt5KPZONv7OPJjBmivaEFvoMdt+6tSp4sz5Yhdm3Z944olijL68vLzY7an97gqZqQGWQEY2NmzFhB4weWe6q+vmdc60Lyb0QXbFQxo28pBAccw8SeYhcwqEAa5oQ2yh37RpE02ePJlOPfVU+s1vfiP+feaZZ8RxtZiIh13xRo8eLc4Yxy56Dz74IE2cOJEn45lRjo2xwsaGLazQqyH1JL15kKouscNnFnpjinogQ2ysD4FezLKbWOiJ6KmnnhJi\/tJLL4m97DHRDnvbNzQ0iPX0L774IuH42m9961t05ZVX0k477cSz7i0r6GmZa2PDFlbok56Ep3LFQp9WyU0mHxvrQzJIZJsqC704CvMdeuyxx8SmOdgRz0\/osTXuJ598QqtWraIpU6YI1rz74mdLpTvhmaxxjJO\/jQ1bGKFPcoMcP9zVTgV+x3a4GKcvdtnIQ7F3svF35sEM1ljot\/GA0P0Pf\/jD3Bn0CNnj1Dqsk1+\/fj2ddNJJNHToULr11lvF1rgI45s2Sc8VMs2oGtGssLFhCyP0YU+Ri4Zi+1Pq4TYs9HHRTP95G+tD+igln6Mr2hB7jH7r1q308ssvE\/71XsuWLaMDDzyQhg0bRn\/9619F+H7ChAn0la98JXmGQuTgCpkhXtm4W21s2FShx8Y0cnKd6k3LyXjqJjZJj8+DXJxvP3VqO81B87SRB+MKswaDmAcNIGpIwhVtiC30hbB8\/\/33xc9yu9uVK1dSz549qWvXrhooIMImPWPHjqV169aJ9MaPH0+TJk3Kpb148eLccEFNTY3ocPhdrpCpBdSMErGxYcsn9Ko3HWTJXRKQ5ztCt1heNvJQ7J1s\/J15MIM1V7QhEaHfvHkzYXc874VZ97h22WWX2CxiyGDcuHFC2HFQjvw8fPhwIejoBGBNP3bnwyX\/77f1ritkxgY1wwRsbNiCCn0am+R4qfMKfZADbZCGjTxkWGwTy5p5SAzaUAm7og2xhB4T7RCy33HHHemWW24RM+4h4tOmTRP\/f++99wj3HHDAAeLv7bffpgULFtBVV12VyDr62bNnCxIh\/vDmV6xYIWb9Y80+fuvTp4+vV+8KmaFKsGE329iwBRX6tMfnQa1X6IMcUctCb06lsLE+mIOePktc0YZYQv+nP\/2JXnvtNTrzzDNzQo+JeM8995zwoOvq6ujvf\/+7mJ1\/8MEH01lnnUW1tbWJnVmvCr36f9Du\/awWBVfI1Fe800\/JxoYtqNCr9wUdK4\/LQL7teYulayMPxd7Jxt+ZBzNYc0UbYgk9wuNPPvkknXfeeXTbbbfRkCFD6M4776RLLrlEsLRkyRLxHcblIfonnHACzZ07NxGhl+P1N9xwgwjlez14ePirV6\/uMIYvi5IkExMFKyoqxP78+OMrPQTQsG3YsIG6deuWSLQniTe5775OdMEFnUTSvXu3UlNTq\/j\/+ed3ovvvb\/u+pWULDRpUJg6akZ+TsMWbpmobflu1aos4srbYZSMPxd7JxmBvm3YAABtOSURBVN+Zh2xZQ1uEv7Vr11JVVZX1R5jHEnqMi8+bN0+E6W+66SYxsx7j5hgHv\/jii8VpdRDez33uc2KsHBPnsMzuuuuuE0vtdF1yfB4CLyfjRRF6aQ+26MUyQL7SQwBDQM3NzWKiJo40tuF66KHdqaqqbWIphH7ZsrbJp\/gOv+F6++13aMSIHtTYWNbhnqTfT7UNecE22FjsspGHYu9k4+\/MQ7asLVy4UAwzywubwEFfbL1iCf1HH30ktr\/FSXUYj8f4fO\/evemwww4TW95eeumlwjv77LPPxN8\/\/vEPMX6PXfTCCj2EG9EAXDjmVo69+4k87okSusewQq9evdijz6A0o1O4ceNGgT3Kjg2X6jXDW4bXjEv16PFdv35t74PT6p54ou2epK+oHr2NPCSNZRbpMw9ZoN6ep\/ToGxsbRTS6pIUehRGCe\/jhh4vJb5htj7D9XnvtRY8++qjYKOfQQw\/NzbKHKEcVej\/avTPt1Xu8oXqejJdtxSmWu41jkkHG6DEJDtvR4gq6O10xrIL87j2qlifjBUHNnHtsrA\/moKfPEh6jJxIz6tHbwb72EHYsqcN3I0eOpKuvvpp23nlnsQ8+vP7TTjtNnHB3\/fXXi85BWI\/eSx06GUgX4\/\/q2nl5Hy+v01fY00jJxoYtiNBj8h020Mla6Hl5XRqlWF8eNtYHfW9vTkos9Nu4wCS8Y489VkzKwxj8T3\/6UyHiGK9H6B5eNyYznH\/++SKsf++994pOwG677RaLTe9mOTIxNazPG+bEgjjVh21s2IIIPU6Nwx+utGbcIy+vR89Cn2pxjp2ZjfUh9ksbmAAL\/TZS7rjjDjruuOPE4TaYlIfwPbz2c889l44++mjCGfRvvPEG3X777fT1r389F943iVNXyDQJ07C22NiwBRF6HCSD+3AFDZ+Hxc7vfnWTHvzOQq8D1fTSsLE+pIdOejm5og2xJuMBboTqP\/30U3r66adpwIABwlN\/6KGHhMhjYhsu7IiH9fX19fViS1r5fXp0Fc7JFTJNwTOKHTY2bEGEHhPw4F2re+FHwSfsMyz0YREz634b64NZCOqxxhVtiC30gBMzE7GPPdbTQ\/QxXi\/3uZdwv\/7663TiiScStqg17XKFTNNwDWOPjQ0bC30YhvneMAjYWB\/CvJ8t97qiDbGFHmfRz5w5k6ZOnSo2m0EBxRj8rrvuSvvss0+OTyytw4Y1WML2hS98wSieXSHTKFBDGmNjwxZE6CUM8OwxRp\/mtcMO7blx6D5N5OPnZWN9iP\/W5qXgijbEEnpMtJs+fTqdc845YokdLsyGnzFjBmHTGfUAGYT4MU4PsceGOXvssYcxrLpCpjGARjDExoZNFXpVyL1nwQMOjNVjeV2alxT6MMMGNvKQJqZp5cU8pIV04Xxc0YZYQg+IMP6unkYHQX\/zzTfFDHvvCXbYYAe7n+FwmR1UdyNjTl0hM2MYY2VvY8MWRugx8\/6aa2JBFPphFvrQkBnzgI31wRjwNBriijbEFnqNmGaWlCtkZgaghoxtbNhMF3ps1INJeezRayigKSdhY31IGaJUsnNFG1joicgVMlMp+QllYmPDFkbo01xDLylioU+osKaQrI31IQVYUs\/CFW1goWehT73y+GVoY8PGQm9E0XHSCBvrg4tEsNA7xKorZNpMiY0NGwu9zSXObNttrA9mIxrNOle0gT169uij1QDNT9nYsIUR+jR3xVOpwYS8oEvr8JyNPGguikYkxzwYQYMzw7os9Cz0RtQoGxu2fEI\/bVr7\/vYS3DBimyUhNvKQJV5J5c08JIVsuHTZow+Hl9F3u0Km0SAXMc7Ghi2f0N95J9HFF3d8YRZ6m0tn+rbbWB\/SRyn5HF3RBvbo2aNPvrYEyMHGhi2f0Kvf49XDLG8LAFWit9jIQ6KAZJQ485AR8J5sWejN4EGLFa6QqQWMjBKxsWFjoc+osJRAtjbWBxdpcUUb2KNnj96I+mljw5ZP6L1nwbNHb0QRs8oIG+uDVQAHNJaFPiBQNtzmCpk2YJ3PRhsbtnxC7z0iNosDbaKWBRt5iPquJj\/HPJjBjivawB49e\/RG1CgbGzYWeiOKjpNG2FgfXCSChd4wVkHI7Nmzaf78+dSlSxdh3eLFi2nKlCni\/zU1NTRs2DBfq10h0zBKQpljY8OWT+jx4uqZTezRhyoKfDPvZ2BMGXBFG5zw6HFc7rhx40ThkELf1NRE1dXVNGfOHPG9\/L96dK4sTa6QaUztiGCIjUKvjsV7xZyFPkIh4EdyCNhYH1ykzxVtcELo4bk\/\/vjjtHnz5pzQ47sVK1bQrFmzqLy8XHj7OB7Xz6t3hUybK5qNDVtQoc\/iLPqoZcFGHqK+q8nPMQ9msOOKNlgv9PDcFyxYQCeccALdeuutOaGHsOOaNGmS+Nf7WS1GkswJEyZQRUUFde\/eXfzxlR4CaNg2bNhA3bp1Ex0zG67lyzvRwIGdhKkDBrTS0qWtObP79SsTR8TiGjWqle65p\/03k9\/NRh5MxjOqbcxDVOT0PIe2CH9r166lqqoqamhooMrKSj2JZ5CK9UIPAT\/++ONzYi5D914PHh7+6tWrc8LvJ\/TyuzFjxtDo0aMzoKN0s9y6dSs1NzdT165dqXPnzlYA8fzzZTRyZA9ha0XFFlq0aH3O7uOO25fWrGnrBAwd+hHV1jZb8U428mAFsCGNZB5CAqb59oULFwoHUl4s9JoBDpMcPPFnn31WiLd3Ml4Uoa+traVevXqxRx+GBE33trS00MaNGwX2ZWVlmlJNNhmE7gcNarMVY\/RPPLEll6GtHr2NPCTLcjapMw\/Z4C5zlR59Y2Mj1dXVsUefFh0Q7rlz54rsBg8eTFOnTqUf\/\/jHBO8bE+z8hB73hgnd295rS4uLJPKxcUyy0Bj9N79JhN9x8Rh9EiXG7TRtrA8uMsJj9BmzirH5sWPH0rp16zpY0rNnT6qvr6cXX3yxQ6ieJ+NlTFiR7G1s2FjozS5TNltnY32wGe98trPQG8aq16Pn5XWGEVTCQn\/ttUTXXGMHHywwZvDEPJjBAwu9GTzkrOANcwwjJKQ5NjZsQT16FvqQhYFvJxvrg4u0sdA7xKorZNpMiY0Nmyr03nF4dYyehd7mkpmN7TbWh2yQSjZXV7TB+uV1Omh2hUwdWGSVho0NGwt9VqXF\/XxtrA8usuKKNrDQ86E2RtRPGxs2Fnojio6TRthYH1wkgoXeIVZdIdNmSmxs2FjobS5xZttuY30wG9Fo1rmiDezRs0cfrQZofsrGho2FXnMh4ORyCNhYH1ykj4XeIVZdIdNmSmxs2FjobS5xZttuY30wG9Fo1rmiDezRs0cfrQZofsrGhq2Q0I8dS4Tz6nHxrHvNhaUEkrOxPrhICwu9Q6y6QqbNlNjYsLHQ21zizLbdxvpgNqLRrHNFG9ijZ48+Wg3Q\/JSNDRsLveZCwMnxGL1hZYCF3jBC4pjjCplxMMj6WRb6rBloy99GHsxATq8VzINePKOm5oo2sEfPHn3UOqD1ORsbNvbotRYBTkxBwMb64CKBLPQOseoKmTZTYmPDFlTo6+vbjqq14bKRBxtwDWsj8xAWsWTud0Ub2KNnjz6ZGhIyVRsbNhb6kCTz7YERsLE+BH45i25kobeIrGKmukJmsfc0+XcbGzYWepNLlN222Vgf7Ebc33pXtIE9evbojaifNjZsLPRGFB0njbCxPrhIBAu9Q6y6QqbNlNjYsLHQ21zizLbdxvpgNqLRrHNFG9ijZ48+Wg3Q\/JSNDRsLveZCwMnlELCxPrhIHwu9IawuXryYpkyZIqwZPHgwzZo1i8rLy8Vn9beamhoaNmyYr9WukGkIJZHMsLFhY6GPRDU\/FAABG+tDgNey7hZXtMFqjx4kXH755VRfX099+\/al2bNni4I0adIkampqourqapozZ474Tv4f93kvV8i0rhYpBtvYsLHQ21zizLbdxvpgNqLRrHNFG6wWegh7nz59fD11ePMrVqzIefiF7nWFzGhF2YynbGzYCgn9tGlth9ng4nX0ZpQxm6ywsT7YhG9QW13RBmuFftOmTTRu3DjhvVdWVm7Hm+rd40fvZ\/UBSeaECROooqKCunfvLv74Sg8BNGwbNmygbt265YZe0ss9Wk7Ll3eigQM7iYdHjWqle+5pzSU0fXonmjGj7bd581rp7LPbf4uWWzpP2chDOsikmwvzkC7e3tzQFuFv7dq1VFVVRQ0NDb46k62VwXO3WugnTpxIZ555Js2cOZPWrVvXYYze68HDw1+9erXoGHgvKfTy+zFjxtDo0aODo8h3xkZg69at1NzcTF27dqXOnTvHTi+NBJ5\/voxGjuwhsho69COqrW3OZVtXtzfV1e0jPuN7\/G7DZSMPNuAa1kbmISxieu9fuHAhLViwIJcoC71efAOnJj36\/fbbT4TncU2ePJl69uwpxDyK0NfW1lKvXr3Yow\/Mgr4bW1paaOPGjQL7srIyfQknmBJC94MGtdnqikdvIw8JUpxZ0sxDZtCLjKVH39jY+P8d9jr26NOiA8I9d+5ckR1m12Ny3SWXXNIhdA\/PHPfNnz+f5s2bJ+6VHnyQ0L3tvba0uEgiHxvHJNUxeozHX3NNOzLqGP0zzxB94xtJoKY\/TRt50I9C9ikyD9lzAAt4jD5jHtDjhQc\/YsSI3NgJSFm0aJHw8B999NEOoXqejJcxYUWyt7FhY6E3u0zZbJ2N9cFmvPPZzkJvAKvqzHqYA+Hv37+\/mIXPy+sMICiECTY2bCz0IQjmW0MhYGN9CPWCltzMQm8IUeqmOOPHj+8w2Y43zDGEpABm2NiwsdAHIJZviYSAjfUh0osa\/hALveEEhTHPFTLDvLNp99rYsLHQm1aK3LHHxvrgDvrtb+KKNli7vE5noXKFTJ2YpJ2WjQ0bC33apaR08rOxPrjIjivawELv0MxKmyuajQ0bC73NJc5s222sD2YjGs06FvpouBn5lCtkGgluQKNsbNgKCf1PfkI0dmzby\/PyuoCFgG\/LIWBjfXCRPle0gT169uiNqJ82NmyFhB6g7rBDG7SffWYExIGMsJGHQC9m2U3MgxmEsdCbwYMWK1whUwsYGSViY8O2ejXR\/vu3AebdMCcjGGNnayMPsV\/awASYBzNIcUUb2KNnj96IGmVjw8ZCb0TRcdIIG+uDi0Sw0DvEqitk2kyJjQ0bC73NJc5s222sD2YjGs06V7SBPXr26KPVAM1P2diwsdBrLgScXA4BG+uDi\/Sx0DvEqitk2kyJrQ2bjRPuCpUTW3mwuez72c48mMGoK9rAHj179EbUKG7YjKCBmAfmwQwEzLCChd4MHrRY4QqZWsDIKBEWmIyA92TLPDAPZiBghhWuaAN79OzRG1GjWGCMoIE9ejNoYB4M4YGF3hAidJjhCpk6sMgqDRb6rJDvmC\/zwDyYgYAZVriiDezRs0dvRI1igTGCBvYkzaCBeTCEBxZ6Q4jQYYYrZOrAIqs0WOizQp49ejOQZx5M5MEVbWCPnj16I+oXC70RNLAnaQYNzIMhPLDQG0LE7Nmzae7cucKa8ePH06RJk3KWLV68mKZMmSI+19TU0LBhw3ytdoVMQyiJZAYLfSTYtD\/EPGiHNFKCzEMk2LQ\/5Io2WO3RQ8hXrFhBs2bNopaWFho3bhwNHz5cCHpTUxNVV1fTnDlzBPny\/3379t2uMLhCpvZSnmKCq1evpgULFtCYMWOoT58+KebMWakIMA9mlAfmwQweXNEGq4Ue3jwu6cWrn9VOQHl5OeE3CIifV+8KmWZUjWhWMAfRcNP9FPOgG9Fo6TEP0XDT\/ZQrPFgt9H4ePUS\/srJSCHu+ToC3MEgyJ0yYQBUVFbrLCqcXAIG1a9dSVVUVMQcBwErwFuYhQXBDJM08hAArwVslDw0NDUJXbL2sFnqALkW6Z8+eVF9fTzI07\/Xg0SlAOEwdw5ekrVmzRohMY2OjrTyy3YwAI8AIMAIJIADnr7a2lnr37p1A6ukkabXQQ8zXrVsnxuhxTZ48mfr37y\/C82GEHs9C7PHHFyPACDACjAAjIBGAwNss8ngPa4RenV0\/ePBgMbnukksuER66DKnAu8d98+fPp3nz5gme\/MbvuQgzAowAI8AIMAKlgoA1Qu8lZNOmTWKWfT6hX7p0aYdQfaHJeKVCNr8nI8AIMAKMQOkhYK3Qgyq\/0D3G6iH+YZbXlR7t\/MaMACPACDACpYKA1UKPtfMYl1+yZIngK+qGOaVCNr8nI8AIMAKMQOkhYLXQlx5d\/MaMACPACDACjEA4BFjow+HFdzMCjAAjwAgwAlYhUPJCr87mt31TBJNLntzvQNrot++BPLPAywPmW4wdO1YspcSKCyynxG6HfAVHQA5zjRgxosPGH4XOg5ATXleuXEmHHnqoWM3SpUuXXKZBz5IIbqX7d\/rx4B2C9A5DMg\/6yoXalnhxxmdX60NJC726HO\/NN9\/MLc1TGzN9Ray0Uyq0YVEhHmQjiP0RhgwZ0mGvhNJGNPjbq0KidqKKTVhVd5f07jRZ7Nng1pXOnfl4gJBPnDiRrrjiityGXyoqzIOeMuJdqSU\/Bz0fxWYeSlroVeLyeTx6ihinUmh5YyEeVEHBrofoFCxatIi9+oBFSnowRx55JL333nsdlqMWOg\/C2yginZkzZ9KNN94ovPowZ0kENNXp2wrx4MVWBYJ5SLZYBD0fxXYeSlboVU8RO+l5PydbvEor9ULYFuNB9fYhMN7PpYVk+LfFXt24MNTh3Xei0HkQ3g6W93OYsyTCW+3eE4V4KFSmmYdky0IhL139zXYeSl7o1TFL3lQnmUqljjHKHGpqasRWxX6RFJUHrwdfyPtJxno3UvXbYKrQNtFenL3h5bBbTLuBYvy38ONBHRdGDup8COYhPub5UpBRlhtuuCF3EJp6wqk63Gg7Dyz0yuQkFvpkKpW3Qqmf0ahhL4R8HS4Wej2csNDrwTFuKvl4kGd2yCO15Wecv6EOmXCHKy4Dbc9LHrB9urpNOgu9HnyNSaVYyNgYQx01RIbFLr300g4T7Ly8cOheTwHIJzBI3e88CNtDlXpQ05+KHw\/eXApNdOQhlPic+Ik8UnV5KKtkPXpJrOzB8WS8+BUoTApq9ET9v5cHb8iMJ+OFQbn93nwhY\/XoZpUHr+foNxkv37PRLCyNp4IKvfTigYo6I595iFdOvDPt1dS8K4Ncqg8lLfS8vC5epQn6tJ9Xfvnll1N9fb1YTsTL64IiGf0+P4EptkTO5uVE0ZFK9kkvD\/kii\/LMDq+nycsco\/MjsVaxVVNzuT6UtNDLSpRvo5boRYqf9CLg3TDHuylOoY2LeMOc+OUpnyfp6gYh8RFLJgU\/Hrwb5ng3heINc\/Rw4d0sR6aq4u1qfSh5oddThDgVRoARYAQYAUbATARY6M3kha1iBBgBRoARYAS0IMBCrwVGToQRYAQYAUaAETATARZ6M3lhqxgBRoARYAQYAS0IsNBrgZETYQSSR6C1tZVee+016tWrF33+85+PlOFnn31GGzdupM2bN4sVDzvssEOkdPghRoARsAcBFnp7uGJLSxwB7x4DGzZsoF\/84he0devW7ZD59re\/Tdgj4sMPPySIu3r99re\/pWnTptHtt99OX\/7ylzv81qlTJ9pjjz24A1DiZY1f3y0EWOjd4pPfxkEE5CYp06dPpzlz5tDOO+8svPqjjjqKXn75ZSH02Lxm1apVdNJJJ4nfIfR77723OMgG58kHvbxLu+Rz6tKk8ePH0wUXXLDdITny3jAHRMlllVjbLPdVCGor38cIMALBEGChD4YT38UIZIaAV+i\/9rWv0cEHH0yHHHKICOHDY6+trRUn1F1yySU5b\/zTTz8VYXp0BPbcc09hPz7j3PMLL7yQvvSlL1FzczMddNBB4p5\/\/etf1L17d1+P3rv1aqEd3sIIPWzypp0Z0JwxI+AoAiz0jhLLr+UOAl6h9x4AJDd8km8sTz+DR79o0SK64447qKqqik444QS69tpr6eOPP6bTTz+dHn74YYInDe8ckYK\/\/e1vosPgN\/7PQu9OeeI3KT0EWOhLj3N+Y4sQyLebF14BIXRcRx55JPXr148++eQTevPNN8XY+\/z586lLly7C229sbKT77ruPzjzzTHFwB75766236IwzzhAH2uC3d999l6688kraZ599fNGJKvQDBw70HT5QQ\/Xs0VtUINlUKxFgobeSNja6VBDA7PjXX3+dMON+7dq1dPfdd4uwO8bou3btKrzy448\/Xgg1xunxf4i5FHrgBGHHb\/Dk8dyNN95ICK\/j+YqKCjGWj9D+fvvtlxfWqEI\/bNiwDmn67TfOQl8qpZnfMysEWOizQp7zZQQCIoDx8JtvvpmOO+444a0jzI4Z97vvvjs99dRTeYUe4v2HP\/yBbrvtNlq\/fj2NHDlSHCD0n\/\/8h4499lh69dVXxd+pp54qZu9jSACT9zDW773yCX2hiX41NTXkFXp0QtSz15EPC33AgsC3MQIREWChjwgcP8YIpIXAs88+S\/fffz9ddNFFdN1119E111xDd911F5111ln0wgsv5IR+zZo1dMwxx+Q8eqyRv\/rqq8Xvp512GmFZHSIChx9+uAjRYx39ihUr6MADDxSdCNyL2foI6e+0004dXk+HR48DQx544IEO0QYW+rRKEedTygiw0Jcy+\/zuxiPw0Ucf0ZQpU+jkk08Ws+yrq6tp6NCh9Mgjj9Ctt95KCxYsEEL+3HPPUe\/evemLX\/zidqF770uqxwLD68faeVyYpY\/Ogd8mOnGFHnmqRxOrNrFHb3wxZAMtR4CF3nIC2Xy3EYD4vvTSS2I5HZbCYcLcK6+8IkRz9OjRIozfv39\/euyxx4TXDtGWY\/SYmIdwfZhLztjHRL5CYhxmeZ2cUHjDDTdQZWXlduaw0IdhiO9lBMIjwEIfHjN+ghHIBAFMzLv++uuF542d7TBGjwvj7PD0scsdwvnYSAdeOSbfISLgvV588UWqq6ujW265hfbaa68OP+M5fLfjjjtqEXo56x4Cjxn+fhcLfSbFiTMtIQRY6EuIbH5VOxHAjPtly5YJkT\/iiCPoqquuEmPsmE2PnfEg+phIh41v4PGPGjVKTIKTIXnvW6uhe6\/nng+hqKF7pIehB7+roaFBePgs9HaWS7baHgRY6O3hii0tQQSwYx2E\/Pe\/\/z1ddtlldMopp1Dnzp3pgw8+oKlTpxK884kTJ9KQIUOEsL\/\/\/vsirH\/22WcTtrP1u3QIvU4qWOh1oslpMQLbI8BCz6WCETAcAYg9ZsF7PfR83yMCgNC7N\/wuXxMhfWx36xeiD+rR64SMhV4nmpwWI8BCz2WAEWAEIiDgPdQm33h72KT5UJuwiPH9jEB4BP4PX5SJxNoodrcAAAAASUVORK5CYII=","height":298,"width":495}}
%---
