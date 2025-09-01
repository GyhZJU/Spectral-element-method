% 参考文献：Spectral Approach in Vibrations of Overhead Transmission Lines
% 目的：对比SEM和FEM结果（workbench项目文件：111）
% 状态：成功
% 备注：比1.0程序更合理。但是在L=33m位置激励时，FRF会丢失第三阶振型。可适当增加n_elements数值来避免振型丢失。
% 原文献图3结果有问题，图中标注张拉力为T=10kN，但其FRF显示的频率远低于T=10kN的频率。
% 猜测图中应该是T=100N的工况。material.T = 100N时可基本复现原文先图3.
% ================================================
%%

clear; clc;

% 材料参数
material.E = 74e9;       % 索的弹性模量 (Pa)
material.ita = 0.01;       % 阻尼
material.rho = 2700;       % 索的密度 (kg/m³)
material.A = 50e-6;       % 索的截面积 (m2)
material.m = material.rho * material.A;        % 索的每延米质量 (kg/m)
material.D = sqrt(4 * material.A / pi());    % 索的直径 (m)
material.L = 100;        % 索的长度 (m)
material.T = 10000;         % 索的张力 (N)
material.I = pi * material.D^4 / 64;  % 索的惯性矩 (m⁴)
material.EI = material.E * material.I;

% 单元划分
n_elements = 3;  % 3个单元，适当增加以避免振型丢失
element_length = material.L / n_elements;  % 每个单元长度

% 频率范围设置
omega_values = linspace(0.1, 50, 100);
fre = omega_values ./ (2*pi());
H_values = zeros(size(omega_values));

% 激励和响应位置设置
excitation_dof = 3;  % 激励节点
response_dof = 3;    % 响应节点

for i = 1:length(omega_values)
    omega = omega_values(i);
    
    % 初始化全局矩阵 (4节点，每个节点有2个自由度)
    global_dof = 2*(n_elements+1);
    K_global = zeros(global_dof, global_dof);
    
    % 组装每个单元的矩阵
    for elem = 1:n_elements
        % 获取单元局部矩阵
        S_local = spetral_element_matrix(omega, material, element_length);
        
        % 计算全局自由度索引
        start_dof = 2*(elem-1)+1;
        end_dof = start_dof + 3;
        
        % 组装到全局矩阵
        K_global(start_dof:end_dof, start_dof:end_dof) = ...
            K_global(start_dof:end_dof, start_dof:end_dof) + S_local;
    end
    
    % 添加阻尼
    % K_global = K_global * (1 + 1i*material.ita);
    
    % 边界条件 - 简支梁
    % 第一个节点(位移=0)和最后一个节点(位移=0)
    fixed_dofs = [1, 2*(n_elements+1)-1];
    free_dofs = setdiff(1:global_dof, fixed_dofs);
    K_reduced = K_global(free_dofs, free_dofs);
    % 设置位移激励值 (单位位移)
    % U = zeros(global_dof, 1);
    % U(excitation_dof) = 1;  % 单位位移激励
    
    F = zeros(global_dof, 1);
    F(excitation_dof) = 1; % 单位激励力
    F_reduced = F(free_dofs);

    d = zeros(global_dof, 1);

    d(free_dofs) = K_reduced \ F_reduced;

    U_in = d(excitation_dof );
    U_ou = d(response_dof);
    H_values(i) = 20 * log10(abs(U_ou));

end

% 理论解计算(前5阶固有频率)
n_modes = 6;             % 计算模态数
omega_b = zeros(1,n_modes);
for n = 1:n_modes
    omega_b(n) = (pi()^2 / material.L^2) * sqrt(material.EI/material.m) * ...
                 sqrt(n^4 + n^2 * material.T * material.L^2 / (pi()^2 * material.EI));
end
f_theory = omega_b/(2*pi()); % 转换为Hz


% 绘制频响函数
plot_FRF(fre, H_values, f_theory); %[output:92ae49f9]

%%


% 谱元矩阵函数
function S = spetral_element_matrix(omega, material, L)
    % 计算波数 k (公式5)
    alpha = material.T / (2 * material.EI);
    beta = material.m / material.EI;
    % 计算波数
    k = sqrt(sqrt(alpha^2 + beta * omega^2) - alpha);
    
    % 计算矩阵元素
    [s11, s12, s13, s14, s22, s23, s24] = calculate_matrix_elements(k, L);
   
    % 计算Delta
    delta = cos(k*L) .* cosh(k*L) - 1;
    
    % 构建对称矩阵
    S = material.EI ./ delta .* [ ...
        s11,  s12,  s13,  s14; ...
        s12,  s22,  s23,  s24; ...
        s13,  s23,  s11, -s12; ...  % s33 = s11, s34 = -s12
        s14,  s24, -s12,  s22  ...  % s44 = s22
    ];
end

function [s11, s12, s13, s14, s22, s23, s24] = calculate_matrix_elements(k, L)
    s11 = -k.^3 .* (cos(k*L).*sinh(k*L) + sin(k*L).*cosh(k*L));
    s12 = -k.^2 .* sin(k*L) .* sinh(k*L);
    s13 = k.^3 .* (sin(k*L) + sinh(k*L));
    s14 = k.^2 .* (cos(k*L) - cosh(k*L));
    s22 = k .* (cos(k*L).*sinh(k*L) - sin(k*L).*cosh(k*L));
    s23 = k.^2 .* (cosh(k*L) - cos(k*L));
    s24 = k .* (sin(k*L) - sinh(k*L));
end

function plot_FRF(fre, H_values,f_theory)
    figure;
    plot(fre, H_values, 'b-', 'LineWidth', 2);

    xlabel('Frequency (Hz)', 'FontSize', 12);
    ylabel('Response (dB)', 'FontSize', 12);
    title(sprintf('FRF'), ...
        'FontSize', 14);
    grid on;
    xlim([0 max(fre)]);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":33.7}
%---
%[output:92ae49f9]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAXQAAADgCAYAAAAT452yAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXW9sVtd5f9is2dY6aWGz\/DeIqbOyfqhoVVW2CGqyqPmQdqiJImpDEI6HHBSpiTWBscmisGiACaaqHFJlhCAHygz5QD\/AlGoCZUuCnFhTkfiULWkTN7XBlpe2IpGwNW9en+uc1+c9vveeP\/fce57Xfq70Cl6\/5977O7\/nOb\/73Of8W7e4uLgIfDADFcrACy+8AJcvX4aRkRFobW1NrcWHH34I3d3dsHXrVujv748t+\/Of\/xwef\/xx2LNnD\/zgBz+Iytjco0JpZNirhIF1LOirxJJrtBo2Ypsm6AsLC3Djxg04dOgQ3L59G1555RX48pe\/zIK+Rv2qUqvNgl6plmPcJbE9efJkIhuDg4PQ0dER\/S4E\/ebNm4nl29vb4cCBA\/DVr361VAYfGqb3YLMwAyEZYEEPyT7fOzMDKLbnzp2DRx55BO66664V19uyZQt885vfLBP0r3\/96yWR\/\/zzz2F0dDSKzo8cOQLf+c53YN26dWXXsblH5grxBZiBDAywoGcgj08Nz4CPlAuK+sGDB+Hq1avw4x\/\/GPAhIB829wjPCCNYywywoK9l66+CutuIbVoO\/Ze\/\/CU88cQTESNy\/hy\/29xjFVDKVahgBljQK9h4DN1ObHWjXH7605\/CwMAAbNu2DZ577jmorq6OKGZBZ0+rFAZY0CvFUowzlgEbsdUJ+meffQZ79+6Fd999F15++eVS6sXmHmwmZiAkAyzoIdnne2dmwEZsdYKOYK5duwZPPvkkYMfp8PBw1NFqc4\/MFeILMAMZGGBBz0AenxqeARuxNRF0HI+Oo11ee+01wCGP3\/\/+9+HYsWPGk5fCM8II1jIDLOhr2fpcd2aAGVhVDLCgrypzcmWYAWZgLTOQu6C\/9957cP78eTh69CjU1taWuJZn3+HEDpyhJ47XX389mq2HhzzTby0biuvODDADzICOgVwFHcV8x44d0WJIsqDj31HQT58+DR988EHp\/+vXr4+mZ+\/fvz\/KW+Ih\/q9beElXUf6dGWAGmIHVzkBugi46q+6\/\/37A4WCyoONveOCKd3fu3InG\/m7fvj2K0jE6HxsbK5XHshs3bixN1V7tBuH6MQPMADPgykBugo5ReJxACwHfvHlzJNLqd1nssVLqd9eK8nnMADPADKx2BnITdDkfLkfcakQuRFtE4WpEjhH7xMRE4vrVk5OTgB8+mAFmgBkomoGWlhbAD5WjogUdhbyvrw\/Gx8ep8Mk4mAFmYA0x0NbWBkNDQ2REPZig+0i5iE5XJLS5uTm4G+GDBWcXMp5kUzBH6W7K\/OibMRWOBA51lJ6+BvmVKFzQ1RRLXKeonGJJ6xQVgk6FUMajd1TmKJ0j5qdyfIiarZC5IILua9giNUIZT+U0RoGUbVZZDxhES8VmVHDIFgwi6CJKF9t6uU4sokYovlmcOXMGurq6oqGWoQ9qeJAPapgYT7qXUuOHkg9R059CIvQ8RY0aoXNzc3Dr1i1obGyEmpqaPKtudG1qeBA0NUyMJ92VqPFDyYeo6Q8LupEsmhei5vzU8FBqjMKq1DhiPPr2RoUjFnS9raxKUCOUiqNRFSsWdL17sw9VDkfU9IcjdL3vWJXgxqinizmqrBQHNXtRCgpY0PXt3aoENUKpOT81PJQaI9W3GGo2o4aHkg9R0x+O0K0eH\/rC1JyfGh5KjZEFXe\/PFO1FCRMLupkPGZeiRig1AaWGh1JjZEE3a2bsQ8k8UdMfjtDNfNq4FDXnp4aHBV3vStRsRg0PJR9iQdf7s1UJaoRSc35qeCg1Ro7QzZoa+xBH6Gae4qEUC3o6idwY9U5GjSPGUzk2o6Y\/nHLR+45VCW6MerqYo8p6CFOzF6W3PBZ0fXu3KkGNUGrOTw0PpcbIKRezpsY+xCkXM0\/xUIoFvbKiPRZ0vdNTE1BqeCj5EDX94ZSLvn1ZlXB1\/l27AH7yE4A33gB46CGrW6YWdsXjD8HKK1HDxHg4KHD1dxZ0V+YSzqNGqKs4rFu3VEFccffjj\/2R5IrHHwIWdFsuqdmMGh6O0NM9Kvf10G0d2qb8ahD0iQmAv\/iL5VovLtowwNFVVraoCRbj0VuUCkfU9CdoygW3lhMbXOzZswf6+\/tLlnz99dfhwIED0ffBwUHo6OiItTI1Ql0c7fnnAf7hH5arhxG6r70xXPDom1O2EtQwMR4OClw9mpr+BBN0FOyxsTE4evQo4J6iu3fvhs7Ozki4P\/zwQ9i\/fz8cO3Ys4ln8v7W1dQXv1Ah1EYe\/\/muAf\/\/35aqNjAA8\/riri5Wf54LHz52Tr0INE+NhQXf1eWr6E0zQMTrHQ0Tl8ndZ7Gtra6GSNol2EQdMt2DaRRws6K7Ny+08F5u53cnsLMaj54kKRyzoX9gqLkJHcW9vb48EPEnsVVNTI9TW0dT8OdYPo3MUdR+HLR4f99RdgxomxsMRus5nk36npj\/BInS8sSCjqakJRkZGQKRU1IgcxR83qpVz7IJgcY3e3l5oa2uL\/tzQ0BB9QhwoDtPT01BfXw\/4dqE7zp2rgp6eqrJijz22AK++uqA71eh3WzxGF81YiBomxqMXdBufzugeRqeHtBlygR88pqamoK+vD9RN7o0qkVOhIKNcULRv3rwZ5dDxGBgYgM2bN0c5dBdBl7np6uqCXTiwO8AxPz8Ps7OzUFdXB9XV1VoE27c3wvj4ys2kP\/rIz9hFWzxawB4KUMPEeNKNSo0fRBsS09mzZ+HMmTNlpK1pQf\/Nb34TdYKKFIuI1lHIT58+DadOnYrIisuvq64nIvShoSFobm6Ofg4ZoWMH78zMTIShpmalUKv4H3ywCq5dq4pGtdx77wL88z8vRet37sx5kE68jh0eLzfVXIQaJsaTbjBq\/Cy1j3B+LUfo4+PjMDw8vLYjdJ2gX7lypSzFslo7ReX8+f33A+BHDF\/0NXSRWn4YGyM1TIwnXdCp8UPJhziH\/oXvxKVcMJeOUflaGbaIQxVxyCIe\/\/ZvSyNduruXvvsa6cKNUf\/OQY0jxlM5NmNB\/8JW+MqEefPLly9Hf1mLE4vkCUUo6HgIgcdI\/eBBvWPrSlATB0rRleCOGkeMR+fVdN7yWND1trIqQY1Qm8YoJhTJ67eINV18DV20wWNFfIbC1DAxHk65uLozNf3BegQZ5eJKoHoeNUJNxUFOt8jiLSYZYT5dRO1ZuDLFk+UetudSw8R4WNBtfViUp6Y\/LOiulkw4z1QcMFf+2mtLF0HhRgHHQ5416mORLlM8nmlIvRw1TIyHBd3V\/1nQXZlLOI8aoSbiII9uUZfLlYXex0gXEzyeTaK9HDVMjIcFXeu0FaI\/HKG7WjJDhC6nW9TOT7mjlAXds3Ey2KwYJEt34QeMnm0qHFELKFnQ9b5jVcLE0eTVFeV0C94I0zA+hy6a4LGqoIfC1DAxHo7QXd2aBd2VuQp55dGJQ1q6BauYFr27UKfD43LNrOdQw8R4WNBdfZoF3ZW5VSLocgSeNNbc59BFamLFKQV9A6BmM2p4KPkQeUEXAJPcjtIiNIiRGqFpzo\/ROaZbxNrnarpFcO5zpAs3RhZQPQMcobtyRE1\/ohz6p59+uoiLZd24cWPFjE25ovLszk2bNkULaa1fv96VCy\/nUSM0TUDlESxPPw0wPBxPgVwu6xIALOh6N6PGEeOpHJtR059I0J9++unFgwcPWokzLrD1\/PPPg+15elPZlaBGaFJjlFMt6lDFuBr7SrtQEwdKr8uCd2ocMR69BlDhiJr+RIK+uOhjCoveCHmUoEZonKOpHaEYdYuJREmcyCNhsgxfLNrxX3kF4IknKusVvmiOdO2A8egYojO0k5r+sKDrfceqhNoY5XHleCHTRbd8DV8sUhxMx9AXicnEeIynsh7AlN7yyAo6bvN24MCByLIiP44bTZw8ebJkbWodogiMGqFCHObnG+H112tK65sj1hdfBHjqKROJWSoj0i5Z1nUpUqzkt4q03H+RmEzYZjws6CZ+EleGmv5EunHhwoXFCxculDo5ca1yFPKtW7dGW8Th3pi4Rnl3dzccP3482siZykGN0F\/8YgFeeukzuHjxT2Bycmn3IcyZ4wJctsvh+lgGoEixEg8grHPaapFFYjLxU8bDgm7iJxUj6NgpKvbzRNBCvJ966qloj09xYBQ\/NjZWEnlXEnyeR0XQMU+O2wxiqkQMSxRijsMTUdRtDznt4rqcblFiJfcTsKDbWrq8fFE2M0VJDQ\/ipoKJiv7Itlz38MMPL8r7e4ot4jo7O8sEHcGLfT99DFeU0zzy2wCCk38bHBwswyGDD0kozup8662l2Z34kQ\/XqFxtZPKYdJfO0aIcX374iDokdbUXhalSBYv50VuOCkch9SeJpSCCjkTs27cPRkZGoLW1NXpQ4EFpCzoRaeO\/+EHxxgP\/rwq4ILelZQH27fsNfPe76+Ev\/3Ip5ZLlyBqlF+X4cv5c1DfpAVQUJlPeGU86U9T44Qg93V5BBD1t42c1tWOySfQ994xCQ0N5br+6GmB+Pr7yclpElJAF3FQMsBxG49hxed99AJ2dc3Dr1i1obGyEmpoam8skls0SpRfVGGWMLOjZzF6UzUxRUsPDgk5M0EVKR07zyBDlaB3\/rn6Xy4pXnl\/\/+m1YWGgx9dFM5TAKv\/fehajjDw95TDk6\/\/T0NNTX10edyT6Oc+eqoKdnKdp\/7LEFePXVBePL5oFHvTl2\/ra2LncAiwfjqVMLsHPnSqxFYDIm6It8rG+b2dxfLcv86NkLyRH6Cn7wmJqagr6+PqA0ArDwCB0Ffe\/evbBt2zY4cuQI3Lx5s2xEjRqRY8Q+MTERpWPUI6ugozjLR3Pz0nfx95aW\/4FHH\/28NGKlvX0u1dvm5+dhdnYW6urqoBpfETwd27c3wvj4UsT\/9tu\/LuHTXT4vPOUP1RrYsaMx+tPo6K3S\/5G3oaHZFRCLwKTjRf6d8aSzRY0fRBsS09mzZ+EMjoCQDnKCjuu4mBw+1nAREfqGDRuiETN4DAwMQFNTUyTaLoI+NDQEzc3N0bUaGhqiT4gD17uZmZmJ7u8r5YL1kKP0LVsW4MoVsyg9Lzwytw8+WAXXri1F6IgLv+OR9DZRBCYb2zOedLao8YNoQ2KSI\/Tx8XEYHh6mFaEXPfU\/LuUij6DBCU14iIjcJOVC5QmZZ75R7ng0nXGaJx4hAwKXWKNGt1pkEZhsBJ3xpLNFjR9ESwUTyVEuRQu6WLVx+\/btpUlKSMz58+ejiP3SpUtlKRaTTtG1IOjqWO+k5Xfl5pm348uYxIxW3YSovDHZiDklcRC4mR+9BalwRFLQcRx6kSkXNJc8kkWkXMTkJpzYtH\/\/fjh27FhkWfF\/HN6oHtQIzdvR5GGMJksC5I0nbocl3To0eWPSy0F5CcbDEbqtz4jy1PQHcZWttigAqpN51HHjrgTI58mTh\/bs2VPW6VkJE4viOChCHGxSL3njidsflQU9W+vI22a26KjhofRWRVrQRSpEdE6qhuep\/\/qmUITz26Re8saj5s8FQ2nrueeNSW8ljtBtOKJmLxb0dOuVInTd+HDfU\/9tnCqpLLUnZFHOL6c60jbMyBNPXP5c2El0jMalhfLE5OJTjIdTLi5+g+dQ05+ylIuI0OWFuuSKYuckjhkXKzC6kuDzPGqEFikOcrojKZ+eJx45taJ20Modo+qaLnlicvEtxsOC7uI35AU9DSCmW06cOFFae8WVAN\/nrWVBRy51+fQ8xSptVyVcWv+LKQagin2emFz8i\/GwoLv4TUUIOoIUqRd55Iu6GqIrAb7PW+uCrsun5ylWIk8el\/KRU0Lq0r95YnLxL8bDgu7iNxUj6K6VC3HeWhd05FwWT\/wuR8R5iZX8IEma5BQ3Agbx5YXJ1f8YDwu6q+9Q05+yHLprpUKeR43QUOIg7+cpR8x54ZHvlzTBKWnp37wwufoh42FBd\/UdavoTCTruWHTw4EGw2bQC0zLPP\/882J7nSlzSedQIDSkOcZ2keeFJGq6o2ikuz54XJlffYjws6K6+Q01\/IkH\/9NNPF3fv3g2YM1cn+MgVFaNgLl++XNpI2uYh4Epa2nnUCA0tDmonaX+\/\/\/XZ04YrqraKi9JDc6RiZDws6K7aRE1\/VqRcBMCkClJZM0Xgo0ZoaHFQO0l\/9rM5uOcevxtuxE33T2sQ8mJdmJ7BJYh9bwLi2iA5p69nLrRPxyGkgoma\/nAOXe\/PViUoOJraSYprlD\/yyF3elvOVx5ibLBAWOdm6JRoxv\/\/++yzoaU5FwYdkfNTwUHoIs6BbyaO+MDVCqTi\/nOpoa5uLRsL4WJ9dfgNIm6Galkvv6JiDwUG\/bw16T0kuQcVmAqEvPA88APDmm1mYWTrXF57sSJavQAUTNf3hCN2nlxFzfpOZpLbVlx8Upmuyi3vIqRffbw229aAcgfoQK9n26kxdW6584LG9p648FUws6DpLWf5OjVAqjoY0YjSNDVvs8WmaHkkzgenolrhrqKkg3N3o299e2t0o5EHJZj4iYrUf5eOPl1Jdrgc1fnxw5MqFeh41\/eEI3Zdlv7gONee\/enV5SziEmKVx24xuSaI1bRKUZ1MYX46azbLikd+ikISREShtaG5MilQwKx6Xe+rOoYKJBV1nKcvfqRFKxdEEjYjnpZc+g76+uuhPJptiJJnApTM07lrydfB3H28Olm5TVpyizbKMApLTLVhRdekFW66o8cMReroFyza4wKLyWi5iDZcXX3wRNm7cCB0dHbb+oC0ftywvb3Chpc2ogGiMu3bdXdrI2TViS1u7xQiM9Bazc+cCXLz4pdJpIUWdmmBlwaOmW7I+xCmJp+xjWTiy8VVdWWoBJeJdsWPRvn37olUVr1+\/DmNjY9FyuTipCCcfdXZ2ehV18fBAIKdPn45mq\/IWdDo3Mv9dOP78fCN85Ss1pRNtUy9yqsT1gSC\/NWAE+sQTd8PVq8s59KzXNWelvCQVcVD5aWxstB6ZpC4BIfpPsnSMUuOH0kOGtKCr66GrOxTh9wsXLpSE17UByefhNd944w24fft26brqfXmTaHem5cb4wgs1gCNTXKK2pIW2XJAlYcJr2Y6ccbm\/eg41wcqCR+60xvQa5tPxsH2AU4yGKWIiLejqjkWqsPresQgj8TNnzsADDzwQrbUuInQUcDz6+\/ujf9XvsmEFob29vdDW1hb91NDQEH1CHNgYp6enob6+Hmpra0NAKLuniufBB6tKqRccZbJly4IW4+RkFbS2LkXSWB7Py3KomFB0nnxy+e0BR2TgPVpast3HFCN1m5nW49q1KkD74vHYYwtRf0lPz9L3U6cWANNcLgc1frAOITFh+8YPHlNTU9DX1weUZtCXUi66CN33jkV4vfvuu68k2rKgy\/l6fLBMTEyUbSItHDNuqYKuri7YtWuXi+9mPmd+fh5mZ2ehrq4OqqurM18v6wVUPO+9VwM7djRGl8UJR+fP39LeAnPdolN1aGgWHn30c+05aQXiOMKHxre+dXfZab29v4Xe3t9lupfJydRtZlIHLDM6+ifw7LN\/HhXHcf74QBScos3Qdi5HKH5kn\/joo4\/LoIfChCDOnj0bBaLyQVLQESAKZFwO\/dKlS3DgwAFvTyK8z1tvvRWJtBr5qykWE0EfGhqC5ubmiOOQETo+FGdmZiIMPmZmujRA+Zw4PA89VBPNHMUD13rBSC7tEFG9mLafBya8JjZgjCgFNvybjzcCHd5KsJmuDvh7nJ1qa5fefDBif\/VVtwg9FD+HDlXB4cPLbxwy\/lCYkEs5Qh8fH4fh4WFvumhiZ12Z1FEu4uSmpiZv28+hMQ4fPgwYSbe2tsYKOt7XJuVC5QmZJf+pM5TL73F45A5O3TBGH2PPVdw6jtRx1PggwaF3Bw+6MKA\/R4dHfwW\/JVzwJNlJnp3r2jHqgicrI+poHdVPQ2CKqxPpHHpWI5iej7nz7u7uaMNp+RAPDRxdI6dYuFPUlNmV5ZIc37ST02QjC1t0Jo0RGzSOV5ejdRR2HN6YZdZjHFYTPLZ1zFLeBU\/SCpjymH\/XjlEXPFnqj+eqD3X8m\/xACoGJBd3QqmrKhYctGhJnUCzJ8dP2\/JQvm2WqfxI8m8b4f\/8H8Id\/WH4l3yNhbPAYUJ65iAuepAe0LIyuw0Jd8GQlQZ0chdeTH0ghMFWkoIuRLmK8uTzBJ6+NonliUVb3Tz4\/zfF1CzjJr73PPgvwj\/\/oB6dtY0yK1jEFg6mYrIctnqz3053vgifpwWv64E7D5IJHV8e03+MmR2F5+YFUNKYkvORTLvJIlsnJySg1gkKO+ey04YNZDJjlXGqEUnE0wWkaHl305nMykWxjV47UdWDwmph+wYau69ilJFg6f7flR9fPIWb4ui4BYItHVz\/d7+oSE\/iwwkPGXzSmihR0dRy6Kpa+x6HrDGvyOwt6Oktpjq+L3vLInyParI1RxiVqj4KOwu6SX8+Kx8RPbcrY4pEfzHFLKGTtGLXFY1PXuLLqEhNx+IvGtCoEHSNyFEwxPjyPmaJZjc+C7i7oeGZaY88jf+5D0PEaGJWisIuZkIIFzK93ddkJOxVxMHmrirN23Ebccjn5AeiSRy+Sn7jO3biO3SIxpbUwavqDWMuGLYq0Sk9PT7R2y4YNG8rWcmlvb4+d4JNVmF3Pp0YoFUczFYe0URBC7G12JjKxo0+O4vLriMGm49QnHpP668rY4jGxk4h68d62wxdt8ejql\/Z73FthXGqwSEwVLejySotiGGFLSwsMDAwAfhdjw7MYzee5LOjZIvSkPLouL5vFhnk0Rozs8OEkFqMS+E6fBvjbv83GUZa6upxrw49sp7Qcufzgto3SbfC41Fc+J+6tMC41WCSmihb0rAYp+nwW9OxiFddpljSu2Yd982yM+IDCKE8Wdt3EpDzxuPBlg0eXP5fv77r8sQ0el\/qKc9KCCNVHi8Kkqw81\/UG8K2aK6ipB6XdqhFJxNGEjEzxxefS8OkQRlwmmrD4W13GaJOxF4LGpjykeeYtBk7SYa5RuisemjnFlZUFXU2aqjxaFSVcnavqzQtDFAl2XL1+OrcumTZu8Lp+rI0z3OzVCqTiajaDH5dF1Y9R1dkn7vUiOkoRdHsNeJB4T3kzxyNG5yZBEWTBNHgA2PmRSL12ZtCBC9dGGhjnIsquTDovp79T0Z4WgqyNbTCsWqhw1Qk0bY1F8meCJy6PnNcKlqAhd5heFDBfHE2vBi99ExN7fT0McbAVUtpHpkgg9PQCvvrp0J906PrZ4svp02mgdNY\/+8ss0bEZNf8oEXZ0lmtVARZxPjVATAS2CF9vGKHKUYjw3vuLaNHqbOoXiKE3Yv\/e938LRo7UkVsg04Uc3hyDNHnL6wmTrPxM8NvZPKqsLImTc778\/B9XVt8BlVycfWMU1qOlPrKDjSBYcnlgJBzVCi3J+U9uY4pGjI4xkRTRrM\/zPNybT69mWSxN2X8sJ2GKSy5vYzHRxtTgc6oxb3aJdJniy1BfPNRlVJb9JdnTMweAgC3oc7ys2uKA4PDHJYVjQ05uSaWOMm1aPVzaJ4Gwbsykm2+valk8T9qzLCdhisRF0E\/HT3V\/dexRFPekowl6mo6rkBxlu4vHII3cFfauipj9lETp+oQgwzTmp4S3C+XWN1UYc5LJxHYi66M0GiyhLjaP\/\/M85+Kd\/uvP7jQruKqtOluUEXHgx5SdLdC7jkq+Tlk8vwl6mo6pk4ccdt\/B7yI1kqOlPmaDLk4qSHJJHufiJiLM0eJtzbRuj3MhtRkLkicnm2i5lBUfz843R3qbyGux4PduJOC4YTB\/CctrBtFMzDY+JqNv6kEv9dfnzpAdRHilBG\/ykBd2mIlTKUiO0COe34d4Wj81uRjY4TAXL9ZpZzlM5ksd3i+vqJidlub96bpLN5GUOfKwyifdVl06Ie0jY+pAtFy4pJHkpgzzSgqZ1oKY\/ZRG6aSUolaNGaN7Ob8u9Kx5sMLZrfphic8Vken3bckl4fK\/qaIorCY\/tuHPT+4Xe7k2ul2nEffXqQrSHqjhCiTo1\/YkV9Lgt4nzuKYo3Ve+xZ8+esnVi5I01BgcHoaOjI9Y\/qRFaKWJl2tjzKFdJHCV1nJoKjwt\/cfyok4JMx52b3j9N1PO2l0ufAGJ66aXPoK+vLqioU9OfFYIuAKoiKgTWx2bM6rrr6vh33oLOtBnqy+XdGPUIVpaghskET+g9TuWZknk9TNQ6ij4UE35c\/ECkfFzmPAhMp0\/fDYcPL0fqeXGTVD\/Sgi6m\/W\/evDk2IkZRHxsbi5bTra2tdbVh7HnybkjqfXiTaHeq82yMrqioYbLBo6ZhfOWyZS7jcvpC9PLqqBb3jxN1nJV5zz35jPnW7ZqV5GMyRy+8UFM2C9hkGQRX31XPIy3oauQcBx7FVWx44YsUvI4s6OpWd2lb3wlCe3t7oa2tLYLU0NAQfUIc6GjT09NQX1\/v\/aHnUh9qeLAO1DDZ4kHR6+mpgmvXliPDv\/u7\/4UjR\/7HxUQrzlHxPPTQ8sibn\/1sLtN2eyYA4+r3k59Mwne\/+8fefVrUraVlAT78cMEEXlRG5Qg78\/Fa4sAH36lTC7Bli\/k1TW+O7Rs\/eExNTUFfXx\/4yFyY3l9XbsXEoqIjdJFPP378eDRDVY3IMWKfmJiIXYtdCLpcya6uLti1a5eu3rn8Pj8\/D7Ozs1BXVwfV1dW53MPmotTwIHZqmFzxDA\/\/adnYdRwXff78LRvzxJaV8fzLv\/xZKU\/s6\/qmALdvb4Tx8WWRRFG\/914\/Dy3EMDlZBd\/61t0RHNu6xdnsvfdqYP\/+uui64njuuU\/h8cdvm1bZqNzZs2fhDC4OJB0kBR3xFZFDl4kQbwXyTkgugj40NATNzc3RpUNG6Ji2mpmZiTCEnPAgOKaGB3FRw5QFD0bp8mgL20gzTkEEnrm5Bti69UuRQOUZcaapGI7Ll7f5w\/QLpjR8HHhdvD4eV67YRdNpNlNTMGiTxx5bWLE4m2sd5Ah9fHz89w\/1YZoRuqhIQl4jAAAT\/0lEQVRgEaNc8F5xYo5\/d0m5UHlC2uRjXR3K5jxqeBA7NUxZ8cTlnbOMQhF4\/vVflyY64VFkXlj1r2efXSjrePQ10cpmMpGKSWezuB2s8phLQDqHbiMUWcumreyopli4U9SdbZ3ju1\/Z\/UxqmHzgSRoh4sIS4vmP\/5iFHTsaS9F5lgeECwb5nLghgllF3WUykYrJZD10m41OXHhiQf\/ilTttj1IetujiWvHn+BArf2iWrkQNk088JlPpdXyqAhoyOpft9V\/\/1VjW8ZhF1E3XbkniysZmea6uWRGCHremy9atW70NV4xL6aDh5HvwxCJdszf73cbxza6YvRQ1TL7xZBX1X\/xiAfAaInceMjpXH8DY8YjYxOE6QzNLusU1KEgSdrwejl\/v6oKor8LmIC\/oCHDfvn0wMjICra2tpbqhwJ44cWLF320qn0dZaoT6FoesnFHD49oYs\/KQdr5vjtT0i20km9cUf1cOVX7UpZZtRT1ruiWrD6UJOwq6zZr41PQHuVkxbHH79u2xG1xgLvvmzZveInVXB5PPo0aob3HIyhE1PFkbY1Y+4s7PiyN5ASnTZYjVjZ9DR+dJ9lLXU7fBmTXd4suHhLDjAxT\/Lx+iA1UXtVPTnzJBx1TL3r174ZlnnimLzkVFEXxeE4tcGyo1QvMSB1d+qOHx1Rhd+ShS0F2Wus2ytZxPTuRrJfmQa2opa7olDx9CQccHjTxEU3CQNjqGmv5YReh5Tv13dUZqhFITUGp48miMrr4jzsuTI1n0TFIvLgtVZa2\/7vw0fmzr5yPdkqcPIb6rV3Em8EpWUNj37gX4m79ZzrVT058yQccvSQB9Ls6lcyCb36kRmqc42PBShFi54MmzMVLEo65imLYksVxWTFCiMDktzadt6of28ZFuKcqHECu+MakbnuD9Ra79r\/7qPdixYwfNiUUmOxbJjYbC7kUs6OkyRu0BU1RjtBH3vDmS0yhpUbpcrrf3t3D0aC2J2cY6fmSR1g2x9JFuKdqH0jpREUtdXR9cuvRobL+jjR\/6KlvqFPV1wSKvw4LOgp7V33SClfX60WvwuqWrpAketQ2Qbd7ycDVI0bGYNOpFjua\/\/W2c7u\/ObBE2i0MXl2uvqpqEd96ZZEF3N+fymSzoLOhZ\/agIcZDFOm7Ei5puefvtX0NjY2NFROjIv0kHsFzHrOuWF2GzNL+SR8j8939fhCtXmmkKujolX57g43NyUdZGKM5nQWdBz+pLRYiDbps1Od3S3z8He\/bks\/64C1em\/HzvewCXLi3dIe6hJdfRduy6itsUk0t9bc6hpj\/R2+Di4nJXjTzWfHJyErq7u6MZnP39\/SsWzbKpeF5lqRFKxdFsXpfzsk3SddcqR\/K4dLVzVI7gcc3zvDaUcLG1qb10m1X46hDFOphicqmvzTnU9KdM0NUNLlSwPA5db2oqjsaCrrdV0RzJ28ipEaoQexw98f77c2Cy8JR5DbOVtPHptL4CXx2iLOjp9ixF6KqgY7SOIi52KML0y4ULF3LZscjV5ag9IW2c37XONudRw0OpMRYt6EmThuTcMnaa4prjlSrocueo+hYixP7++wHwgZbloOLX1PQnNuWCf+zp6YHdu3fDhg0boqn+uKA8fpc3oshiEF\/nUiOUiqMVLVY29lzLHAnBk\/cGVTsLMYdeqYKe9Bbis0OUUlBATX9WCLo8Fr2pqSlajKulpQXSlru1acy+y1IjdC2Llalt1zJHcq5cRLCrSdCTFhYzHYtfaT5ETX9WCLopoVTKUSN0LYuVqU+sZY7kCFYIuix2OJyvkiP0SFBixtz77BDlCN0wh27aIIsox+uh+2GZmnhSaowh0lKyoIuhfWr02tlZuSkX5DSujj47RCn5ELWAMjZCl9MuYuz5iy++CBs3boSOjg4\/SpNyFd6xyB\/FLOh6LovkSI5UhaDLaQrsLGxvr2xBjxu+GNd3oLdMcokibZaGk7ygyxtcXL9+HcbGxso6RTs7O3MXdXVVR95T1N31qTi+XANqmIrEo4o3jvhQ0xGVLuhq2gU3jEBBx8PHCBeO0A1TLjiSBTs\/N2\/eHIm2KqxFDVtEAccDJzPhoX6Xq0PtCVmkOJjIPDU8lBpjiJRLXPS6GgVdHr6I\/QL4wSPrlP8QNqvYCF0dh64KelETi9SIHHFMTEyUBD5O0Ht7e6GtrS36qaGhIfqEOFBAp6enob6+Hmpra0NAKLsnNTxC0NcqR+fOVUFPT1Vko1OnFmDnzgU4dKgKDh9e+htOKmpoqHwfunjxD2Dnzj9a4f84Cxaj9KxHSL9G38UPHlNTU9DX10dz+VxdhF7UFnQugi47SFdXF+zatSurzzidPz8\/D7Ozs1BXVwfV1dVO1\/B5EjU8WDdqmIrEg5ss79jRGJl4aGgWHn30c+jrq4OLF78U\/e2jjz5eNfwMD\/8pDA\/fVebOWD8fR5E2U\/GePXsWzpw5U\/bn0dFRmotzJeXQL126BAcOHCjkSeSSchkaGoLm5ubgETo+FGdmZqI3BAqbE1DDgwaihqlIPDjm\/CtfqYn89O\/\/fgGefXYBHnywCq5dW4rQ79yZW1X8YL2wfniIZQ18CHqRNlPxyhH6+Pj47x9aw4XooilvK9ZDj9voQkwyam1tNb2uczk1xcKdos5UklnESK4Btbx+0XjUcdrqkL6i8ei8KysefIihmPs8smLyhYVaHx7Wy2qDizfffBO+9rWvwfr1631xsuI6PGzRH7VUHJ8FfZkBVdDVIX3UbEYNDzJJBRNZQZcn8sRtLSfy65988kkhi3PxxCI\/ok7F8VnQVwq6GMInBF18p2YzanhY0NO1Yd277767uG\/fvmjdFkypqJ2f4imElxkcHMx9HLqNlFF7QlJzfmp4KDVG4WdFc6QKuLoKYdF4dO2NGh5KPkRNf6KUy9NPP72IOXIx7lsevvirX\/0q6gwtMoeuczD5d2qEUnN+angoNcbQgo55ZZwZKibdiP1GqdmMGh5KPkRNfyJBf\/jhhxflGaAivfK73\/0O3nnnnd9vh7Undgy4jfDmVZYaodScnxoeSo0xlKCr+4sKQReTbqjZjBoeSj5ETX9Kgo7ROa51Lg5Mu5w8eZJcikV9MFAjlJrzU8NDqTGGEnR18SoWdPtwj4pfU9OfVEG\/efNmtI4LhRmPSSanRigVRwslViZNc61zJAv6yMjS6oR4cIRu4j1LZaj4EDX9SRV0\/FHk1c2pLrYkNUKpOBoLurkfFm0zee0WWdDx\/5hHLxqPjilqeFjQ0y0W5dDjUi4s6DpXX\/k7NeenhodSYwz10JMFXV64SmwcTc1m1PBQ8iFqASVH6PaanXoGNeenhodSYwwl6PKKizj2HDe4wIMF3bwxUvFrsoJ+48YNIzbjJh0ZnZhTIWqEUnG0UGJlYua1zpEs6Dh0EafGs6CbeM5yGSo+RE1\/ogh9cVHsbmhHKoXS1Ail4mgs6ObeWbTN5C3nZEEXOxgVjUfHFDU8lN7yqOkPC7rOmy1\/p+b81PBQaoyhHnqyoMvuJcIqajajhoeSD7GgWwqkrjg1Qqk5PzU8lBpjKEHHFIsYe86Crmvh8b9T8Wtq+sMRups\/JZ5FxdFCiZUJncwRgFi\/RfCFqRdMufADz8SD6AztZEE3s5dxKWqEsljpTcccsaDrvSS9BBUfoqY\/HKFn9SzlfCqOxhG6uWFD2EzeRBmRiqVzOUI3s1sIm8UhY0E3s5dxKWqEUnE0FnRjFwoyM5MF3dw+cSWptDNq+hMsQsddibq7uwHXi8FDXdGRN7jI5vAs6Ob8hRAHecVFRCqWzuUI3cxuIWzGEXqCbeT11nGFR\/FdLOHLW9CZObVJKSqOL2OlhikEHnmBLuRGLMzFgm7i1dwpmsYSiYlFuFwvHrimDEbnY2NjpZUeeZNoMyen\/GrKgl5uHRZ0d5+m9NDjlEuCHWVBl\/+PxdXv8iWoERoi2ktrGtTwUGqMIdNSLOgs6NkYSD47eIQu8unHjx+PNtlQI3KM2CcmJmKX8hWC3tvbC21tbVEtGxoaok+IAwV0enoa6uvrSawjTw2PEPS1ztGhQ1Vw+HBVyUVPnVqAnTsXou\/UbEYNT2iO0Hfxg8fU1BT09fXB6Oho2QZBIbRH3DOooIv8OQq5WHvdRdBlAru6umDXrl1BOJ2fn4fZ2Vmoq6uD6urqIBjkm1LDg9ioYQqB5+LFL0FfX13JVKOjt6C9fS76HgJPmqNSwxOao7Nnz8KZM2fKKFtTgi62s0MGtm7dWsqNx4k5lnFJuQwNDUFzc3PwCB33Y52ZmYneEGpqaoILOjU8SAg1TCHwnDtXBT09yxH6lSsLsGXLUoQeAk+ao1LDE5ojOUIfHx+H4eFhjtDVkS2yQ6kpFu4UdX8ucA5dz10IjuQldBGhWAtdpBNu3boFjY2NJIKCEPzorEYFE7U+POSt8JQLPvEHBgagqakpNi\/OwxZ17mz+OxXHlxFTwxQCj7riolg6lwXdzLdD2CwOGQs6AKiTigRRcjqGJxaZObauFBXHZ0Evt5S64qK8IwE1m1HDQ+mhx4KuUyDL36kRSs35qeGh1BiFq4XiSF5xkQXdruGHspmKkpr+BEm52JkuvTQ1Qqk4WmixSrMac7TEjhB0eelcfuCZqQMVH6KmPyzoZv5jXIqKo7GgG5ssyOJcLOjm9okrSaWdsaBns+OKs6kRSsXRWNDNHS2UzcSKi\/LSuRyhm9ktlM045WJmH+dSLOjp1FFxfBklNUwh8WDaRd2iPSQeytEwRR+ipj+ccnF+lMSfyI1RTyhzVFkPYWr2ovQWw4Kub+9WJagRSs35qeGh1BippqWo2YwaHko+RE1\/OEK3enzoC1Nzfmp4KDVGFnS9P1O0FyVMLOhmPmRcihqh1ASUGh5KjZEF3ayZsQ8l80RNfzhCN\/Np41LUnJ8aHhZ0vStRsxk1PJR8iAVd789WJagRSs35qeGh1Bg5QjdrauxDHKGbeYqHUizo6SRyY9Q7GTWOGE\/l2Iya\/nDKRe87ViW4MerpYo4q6yFMzV6U3vJY0PXt3aoENUKpOT81PJQaI6dczJoa+xCnXMw8xUMpFvTKivZY0PVOT01AqeGh5EPU9IdTLvr2ZVWCmvNTw0OpMXKEbuba7EMcoZt5yhel8CmH28ydPn0a1q9fH\/21Uje4mJiYiDaQxY2qN+K6qIEPaniQDmqYGE+6k1Ljh5IPcYSu+I7YWxT\/LAS9krego2ZganjQztQwMZ50QafGDyUfoshN4XuKyu6Dkfgbb7wBt2\/fLgk6\/m1sbAyOHj0KtbW1UfSO0W5HR8cKz6NGKOPRv5IwR5UloNTsxYKe7j\/BBB0jcUxPPPDAA3DixImSoKOA49Hf3x\/9q36Xq0PN2RgPC7qeARb01cIRtfaOvAYTdBTq++67ryTaIuWiRuQYsWMeTwh8nKD39vZCW1tbVj\/JfP7U1BT09fUB40mmkjlKdzPmR98MqXAkcIyOjkJ7e7seeAElggg6PtneeuutSKTVTlEbQZ+cnIwEdHx8vACq+BbMADPADJQzgIHk0NAQtLS0kKAmd0FHgT558mRU2a1bt8Jzzz0HP\/rRj6KRIK2trbGCbppywXIo6vjhgxlgBpiBohlAIaci5kFSLpg77+7uhps3b5Zx39TUBCMjI3D9+vWyFEtap2jRxuP7MQPMADNAmYHcI3Rd5dWUi82wRd21+XdmgBlgBtYSA+QEHck3nVi0lgzFdWUGmAFmQMdAcEHXAeTfmQFmgBlgBswYYEE344lLMQPMADNAnoGKFXR59AylcaB37tyBgYEB2L59e9CxqWrn8549e2LH8hfloWKZhxs3bkS3HBwcjJ39WxQe+T5x6wmFwKFyRIEnOf2Jo9TEDO6i+YnjBjFs2rSpbB2oonFRu19FCrrcAD\/44IMVi3uFIlmI+eXLlyHkQ0Y4P47zxwkP4ntnZ2cQERW8bN68Obq\/eNgcP3486EMP\/SRuPaFQ\/oO8HDlyBH74wx+WFqoLhQXvi+1s37590egzHGKcNmu7aJyqTxV9f6r3q0hBlx2LWkT8jW98Az755JMoGqYyewydjxtjfBOMW08oVGNFAT1\/\/nywKFitN+Uhw+qaT6FsRu2+FSfo6pOZypMapwHjgQuK7d69mwU9xdPloakY+YU6ktYTCoWHkkipb3mhOIm7L2VsoXmqWEGXc9SUIgmKzkYlxSGnpELn9MVbS9x6QqEapdwvhBhC56z37t0L27Zti9JAOBEwJB7ZJpQefKF8Jem+LOieLUJN0AUeTP\/ELXDmufpGlxPCjrODQ2FKW0\/IqBKeC6mchOZI+M2GDRuiFBAe2Nkf0maIgcobuWfze7tcxQq66GCjZmBKgk5RzIXnhhxZgj5z+PDhxPWEvLWujBcKyVGcH4fEI6hEXPjm8Mwzz0QdtXyUM1Bxgi5elcWmF1Q6RWWHo5BDDz2yRdfQQr4269YToiIUITtJ49pVSDwUAgGdT1P4vSIFneqwRTQohQg99Ou66tjqWxS1hw2VyFMOBChwJD90RcpFvBmHEq+0\/RFCYaJ034oUdBGli2V5Q475Vo1JQdCTItCQnVrqxBAKnaLUoj6KHMkTiyjYjNLwW0pCLrBUrKBTJJMxMQPMADMQkgEW9JDs872ZAWaAGfDIAAu6RzL5UswAM8AMhGSABT0k+3xvZoAZYAY8MsCC7pFMvhQzwAwwAyEZYEEPyT7fmxlgBpgBjwywoHskc7VfSl1rRK4vr0ttbn15LLW6RK18Fdv5BDyL0twGq7UkC\/pqtWwO9UJBx0WaQm1ykEOVCr+kuua5T0HHylCYzVk4qXzDEgMs6OwMxgywoBtTFVswbt0h34JObSmMbIzx2bYM\/D+SGhD0mvTKggAAAABJRU5ErkJggg==","height":0,"width":0}}
%---
