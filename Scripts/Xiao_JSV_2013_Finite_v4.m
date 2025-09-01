% 目的：找到FRF求解方法
% 文献：Flexural wave propagation in beams with periodically attached vibration absorbers: Band-gap behavior and band formation mechanisms
% 
% 状态：不需要继续
% 备注：将Xiao文献中的铝梁参数换到Dutkiewicz文献中，观察能否成功计算FRF
% ================================================
%%
%[text] ## 参数设置

clear; clc;
%%\
% 
% %% 1. 几何参数 (Table 1)
% geometry.h = 0.015;      % 梁高度 [m]
% geometry.l = 0.01;       % 间距 [m]
% geometry.r1 = 0.013;     % 内半径 [m]
% geometry.r2 = 0.015;     % 中间半径 [m]
% geometry.r3 = 0.018;     % 外半径 [m]
% geometry.L = 0.1;        % 梁总长度 [m]
% geometry.N = 6;         % 单元个数
% 
% %% 2. 材料参数 (Table 2)
% % 铝层参数
%     material.aluminum.E = 7e10;          % 杨氏模量 [Pa]
% material.aluminum.G = 2.7e10;        % 剪切模量 [Pa]
%     material.aluminum.rho = 2700;        % 密度 [kg/m^3]
% 
% 
% %% 3. 截面特性计算
% % 梁的横截面积
%     geometry.A = pi*(geometry.r3^2-geometry.r2^2);
% % 梁的惯性矩
% geometry.I = pi*(geometry.r3^4-geometry.r2^4)/4;

%%\

% 几何参数
    geometry.r1 = 0.013;     % 内半径 [m]
    geometry.r2 = 0.015;     % 中间半径 [m]
    geometry.r3 = 0.018;     % 外半径 [m]

% 材料参数
    material.E = 7e10;       % 索的弹性模量 (Pa)
% material.ita = 0.01;       % 阻尼
    material.rho = 2700;       % 索的密度 (kg/m³)
    material.A = pi*(geometry.r3^2-geometry.r2^2);       % 索的截面积 (m2)
    material.m = material.rho * material.A;        % 索的每延米质量 (kg/m)
    material.D = sqrt(4 * material.A / pi());    % 索的直径 (m)
    material.L = 0.6;        % 索的长度 (m)
    material.T = 0;         % 索的张力 (N)
    material.I = pi * material.D^4 / 64;  % 索的惯性矩 (m⁴)
    material.EI = material.E * material.I;

% 单元划分
n_elements = 3;  % 3个单元，适当增加以避免振型丢失
element_length = material.L / n_elements;  % 每个单元长度

% 频率范围设置
omega_values = linspace(0.1, 12000, 12000);
fre = omega_values ./ (2*pi());
H_values = zeros(size(omega_values));

% 激励和响应位置设置
excitation_dof = 1;  % 激励节点
response_dof = 2*(n_elements+1)-1;    % 响应节点

% 检查波束随频率变化，后续画图观察
k_record = zeros(size(omega_values));

for i = 1:length(omega_values)
    omega = omega_values(i);
    
    % 初始化全局矩阵 (4节点，每个节点有2个自由度)
    global_dof = 2*(n_elements+1);
    K_global = zeros(global_dof, global_dof);
    
    % 组装每个单元的矩阵
    for elem = 1:n_elements
        % 获取单元局部矩阵
        [S_local, k_record(i)] = spetral_element_matrix(omega, material, element_length);
        
        % 计算全局自由度索引
        start_dof = 2*(elem-1)+1;
        end_dof = start_dof + 3;
        
        % 组装到全局矩阵
        K_global(start_dof:end_dof, start_dof:end_dof) = ...
            K_global(start_dof:end_dof, start_dof:end_dof) + S_local;
    end
    
    % 已知条件：左端点位移 = 1
    known_disp_dof = excitation_dof; % 位移已知的自由度
    unknown_dofs = setdiff(1:global_dof, known_disp_dof); % 其他自由度
    
    % 分区处理系统方程
    K_kk = K_global(known_disp_dof, known_disp_dof);
    K_ku = K_global(known_disp_dof, unknown_dofs);
    K_uk = K_global(unknown_dofs, known_disp_dof);
    K_uu = K_global(unknown_dofs, unknown_dofs);

    % 已知位移条件
    U_k = 1; % 左端点单位位移激励

    % 求解未知位移
    % K_uu * U_u = -K_uk * U_k
    U_u = K_uu \ (-K_uk * U_k);
    
    % 组装完整的位移向量
    U_full = zeros(global_dof, 1);
    U_full(known_disp_dof) = U_k;
    U_full(unknown_dofs) = U_u;

    % 提取输入和输出位移
    U_in = U_full(excitation_dof);
    U_out = U_full(response_dof);

    T_N = abs(U_out / U_in);
    H_values(i) = 20 * log10(abs(T_N));

end

% 绘制频响函数
plot_FRF(fre, H_values); %[output:8f43c4bb]

plot_k(fre,k_record); %[output:22b1053c]
%%


% 谱元矩阵函数
function [S,k] = spetral_element_matrix(omega, material, L)
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

function plot_FRF(fre, H_values)
    figure;
    plot(fre, H_values, 'b-', 'LineWidth', 2);

    xlabel('Frequency (Hz)', 'FontSize', 12);
    ylabel('Response (dB)', 'FontSize', 12);
    title(sprintf('FRF'), ...
        'FontSize', 14);
    grid on;
    xlim([0 max(fre)]);
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
%[output:8f43c4bb]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAUEAAADCCAYAAADXXXDzAAAAAXNSR0IArs4c6QAAIABJREFUeF7tXQ+Q10XZf6gz7qbX9018EcETSYeKsQaLCmSUrLRJDZk0g8MJPPFAE7kKjn9qaCqcHDUdaIpAJ5cCOqPG0UijSCJ0wiRMZuiMlGKeIHNppU4CXvK+z\/fn\/m5\/e98\/u\/vd\/X73t79nZ5jh7p7dffazz36+z\/57tt+xY8eOASVCgBAgBCoUgX5EghXa89RsQoAQCBAgEiRDIAQIgYpGgEiworvfXuPvuOMOWLlyZWwFS5YsgYkTJ8K+ffugvr4eDhw4ECp\/yimnwLnnngvXXXcd4P9ZUqnDXkup5HJHgEiw3HvQUf2RoO6\/\/3749re\/DSeccEKolueccw586UtfKpLg5z\/\/+YAU+XT06FHYvn07rFu3DvDvy5cvh5NOOikQUanDUZhILQcQIBJ0oBN8VAEJatOmTdDW1gbDhw+PbSLzBMePHw\/z5s3rI4t7dw899BAsWLAAmPfISFC2Dh8xpjaZQYBI0AyOVIqAgEkSxKLDiFKlDuogQiAKASJBsg0rCKgQVJIniAru3r0brrzySpgxYwbMnDmzOB0mT9BK91VUoUSCFdXd2TXWFAn29PTAc889B7fddhu8\/fbbcO+998IZZ5xBJJhdV3pfE5Gg912cTwOTdm75tb2k3WFswZgxY4I1wc997nPFBqnUkQ8KVGs5IEAkWA69VIY6Ju3csp1hfr2P3x1+9913gx1h9AIXL14MF110EfTr168ECZU6yhBCUjkjBIgEMwK60qoxMR1GIly0aBFs2bIF7rrrLkDi5JNKHZWGP7VXHgEiQXmsSFIBARWCitsY+etf\/wrTp08PaubXA\/FnlToUVCfRCkOASLDCOjyr5qoQVNLu8COPPALz58+Hyy+\/HH784x9D\/\/79g2ao1JFVu6me8kOASLD8+qwsNFYhqCQSfOedd2D27NnwzDPPwN13312cFqvUURagkZK5IEAkmAvs\/leqQlBJJIho7dixA6699trg6lxra2twFU+lDv8RpxbqIkAkqIsc5YtFQIWgZEgQzwviLvF9990XXJ377ne\/C0uXLpW+mkfdRQhEIUAkSLZBCBACFY0AkWBFdz81nhAgBIgEyQYIAUKgohHIjAR37twJkydPDsDGS\/B8yCR2\/WnIkCFSoZcquseo8YQAIWAUgUxIEBe+586dGyxk19bWBme+xo4dGwTQfPDBB6GzsxOam5uDK1JIiGvWrIEBAwYYbSgVRggQAoRAGAKZkCB6gTy54f8xoTeI\/x82bFhAiG+99RZMmzYt+D1emKdECBAChIBtBDIhQZ7cRo4cWfQEL7nkkuD\/dXV1Aem99957JV6i7cZT+YQAIUAIZEKCCDMjQpzyYnQQJD3R84sjwa6uLsB\/lAgBQsA\/BHCZDP\/lkTIhwajp8KxZs6Q8QSS\/pqYm2LVrVx4YUZ2EACFgGYHRo0dDS0tLLkSYCQnya4CIJU+Kq1atSlwTZDvLCBL\/5KLlfnGqePwA4HUxwoAw8M0OmG2zGWLWAy8TEgzzBPGNWdwR7ujoSNwdZiSYF0hZd0pYffv374e1a9fC1KlTg49GJSbCAMBHDPIe35mQIA5YPAqD4dEx4eYIfwwm6Zxg3iC5QDiHDx+GgwcPwuDBg6G6utoFlTLXgTAA8BGDvMd3ZiSYZsTkDVIa3U3l9dH4VbEhDIgEVW1GRp5IUAYlB2SIAPwkAFXT8tEO8nZyiARVrTAneR+NXxVKwsDPDwGRoMRIyBskCRWtixAB+EkAqobjox3kPb7JE1S1wpzkfTR+VSgJAz8\/BESCEiMhb5AkVLQuQgTgJwGoGo6PdpD3+CZPUNUKc5L30fhVoSQM\/PwQEAlKjIS8QZJQ0boIEYCfBKBqOD7aQd7jmzxBVSvMSd5H41eFkjDw80NAJCgxEvIGSUJF6yJEAPYIoF+\/QvcdPQpw3HHWuzJVBT7aQd7jO9QT5EPhh\/VY1nd48wYpldUayuyj8atCYwOD\/fsBPvnJgiY33wywaJGqVtnK28Ag2xb0rS3v8V0kQT7en\/gGCK82i\/m3adOmPneAbYGZN0i22qVSro\/Gr9J+lLWBAZGgai+Yl897fAckiAR4yy23wKJFi5Te9tDNpwpj3iCp6mtD3gYB2NDTZpk2MCAStNljcmXnPb5pY0Sun3KXskEAuTdKUQEbGBAJKnaCBXEiQQlQ8wZJQkXrIkkEwBb3jx2zrkpuFSRhoKMYkaAOambz5D2+SzzBsJh\/GPl55cqVxVZnvSmCFecNktku1ystjgDuuw+gvr5Q7u9+B3DeeXp1uJ6LSNDOumje\/Z73+C6SIBLghg0bisFOWaDT8ePHBxGga2pqAN8Prq+vh2XLlmX6JGbeIOVtJFg\/kaAdAiBPMH\/rznt8ByQY9sobI7zrr78+eBOYJf6xdCTGLFLeIGXRxqQ6iATtkCDizpYS6IhMkhXa+Xve47u4Oyw+es6OzEyaNKmEBFFh\/iF1WVgYqeLbIq6E18fzYegJ4JMdr7wi25J85IgE7ZAgeYL52DNfa0WQoPi+MHqT+GDMvHnzgrdHOjs7gyk3vkkcRrC2QCISzH8AqGhAa4J2PgQqfWBD1tb4ltU1E08QG7lt27aA9MSEpIevp+GUWyRLJmsLJCJBWTNxQ842CV55JUBbmxttjdLCBgZ5t9jW+JZtVyYkiN5ed3c3PPnkk4G3x6bDuKY4f\/58qKurCzZawtYmsSEMpMbGRpgwYQKcfPLJsu2LlRsxojqYDtfW9sC+fT1GyrRVCBr\/G2+8AYMGDQo2qfh0\/\/1V0NBQFfxq8+bDXu8OR2Ggi3tXVxUMH17A7ooremD16vK1A10M8syH\/bl7925oamqCPE6eYNszI8EVK1ZAW1sbDB8+PJjyYmpoaAB+LTKJBDEPvrs7ZcoUI\/02btypgIMASfDpp18zUqatQo4cORJ8SAYOHAj9+\/cvqebhh\/8LmpoGBr9bt+4gjBlz2JYauZYbh4GuYtj\/aAeYLrvsXWhp6dYtKpN8NjDIRPGIStrb24P3tAu2uy7TUydMpRISRC9NJokbG0l5+DVAlGWbK3feeScsXbpU2hNsaWmBUaNGGfcEcWPkxRfdJg78QBw6dChou\/juMO8JPvFED5xzjtveTJK9RP09DgPdMnEmgDMCTOXgCdrAQBc7E\/nQE9y4cSO0trbmS4ImGhNXBpLe+vXri+cN+R1mPIxNa4LJPUC7w3Y2BfjdYVoTTLZDGxJOrAnaaBhfJm54zJ49GxYuXAi1tbXBOuCQIUNod1gBeCJBIkE0F9oYURg0kqKZTIdRF\/6cIH8LBf\/GbqcgMbJ1Q15\/W18K2h2WtBJHxGwQAHmC+XeurfEt27LYoKpLlizpc1B6zpw5oUQlW6GOnC2QiAR1eiO\/PESC5AnasL4+JMh2aNl0VazUp2tzPpIgnnPDtS0fE5EgkaANu+5DglEHllnlutfm0ihPnqB8AAUiQTVLo+mwGl42pG2Nb1ldIz3BsWPHlkyFWYG4fof3f1lkGdmK0sjZAok8wTS9kn1e8gTJE7RhdbFrguLhRZwK84eebSgUViaRIHmCaBdEgnYwyGocR9Vja3zLtisyvD7\/8BIrTNzVla0krZwtkMgTTNsz2eYnEiQStGFxFf3GCJGgDZOyVyaRIJGgDesiEvQsniBtjKgNE9oYUcPLhrStmZ6srhX95CZ5grJm4oacDU\/wqacAvvrVQvvo2lw+\/ewECWLTK\/HxdSLBfIxet1YiQZoO69pOXL7Y3eGojFmHvLH1pfCRBOm1ObVhQp6gGl42pG2Nb1ldaU3QgzXBW24BwEeCMNGaoKzpF+SIBNXwsiFNJCiBqi2QfPEEiQQljChChEhQHztTOW2Nb1n9yBMkT1DWVnKXozVBWhO0YYREgp6RIK0Jqg0T8gTV8LIhTZ6gBKq2QKLpsAT4DonY8ATpnGD+HWxrfMu2jDxBzzxB2hiRNX3aGFFDyp60syTInxtkd4aXL19efA9EFxKMMD137tzggSV8eQ4TRZZORjPOC6KNkWT8oiRoOqyPnamcTpIgKsUiSO\/Zswc6OzuD0FkYcBWfyJw0aVJomK0kUFjAVnxnlIXR54O04mt3SIhr1qyBAQMGFIuzBRJNh5N6zK2\/03SYNkZsWGRiPEExkjT+vGHDhj5EJaMc5n3++edh7969RU8QSY9em0tGjx5askMAvCd43nkAuLHkcrLxIci7vbacHNl2JUaWFklQN7I0ToMXL14MM2bMCLw9nA6zl+fq6uqCR5eTHl83fVOFPEFZM3FDzgYBEAnm37fOkaBIRCIJ6kaWxnxf+cpX4MQTTyyuCeL\/cXo9b948KRJsbGyECRMmGH98vba2B\/btc\/vBciQAfKh60KBBUFNTU2K5lfL4ehwGukMZSfDCCwuPr+Oj9fh4vcvJBgZ5thdtGpfHmpqa3Hp8PWpNsKOjAxYsWKCsLJa3bdu2gOz4jRFVTxA7a+rUqTBlyhQj\/TZu3KnQ1VUFSIJPP\/2akTJtFXLkyBHo7u6GgQMHQv\/+\/UuqaW39BLS2nhD8rqWlGy677F1bauRabhwGuort3FkNkycPDrKPHn0Y1q8\/qFtUJvlsYJCJ4hGVtLe3w9q1a4O\/mp7pybZLKbJ01LvASZWx3V9ejpX1yCOPSK8JtrS0wKhRo4x7gqjXe+8dTmpGrn9HD\/3QoUNB26urC54LS7fdVgW3314V\/HjDDT1w441uezO6QMZhoFvmjh1VcMEFBexwTXDz5vK1A10M8syHnuDGjRv\/\/yPe6h4J2gJGPCLjwu4wtvXYMVstNlNu3HrY6tUADQ2FeuicoBretCaohpcNaefWBG00ki\/TxXOC5U6CCxcCLFlCJKhjuzwJDhsG8MorOqVkl8fG5lB22ofX5CQJsoPS7Dwgemu4Fogpj8eWbIHUr19vp5SzJ1hfD3DffYW2YEitRYvyNms79dsgACJBO32lUqqt8S2rQ+iaIL8D3NXVBfX19QH54cYG\/g0T\/j+rZAskIsGsetBMPUSCds5Kmukd\/VJsjW9ZjRLPCYoK6p4TlFUoTM4WSESCaXol+7xEgkSCNqwukQTR80MSYlfZ0twY0W0AkWC88fPT4ZtuAvjJT3SRdjufbRLE1pfzsojbvRetna3xLYtH5HQYC2hoaAgOMw8dOrTk7jDe7qDpsCzEZuTiCABfS8O1rcIyBUBzs5k6XSvFBgk+9hjAxRf3tpRIMPted5IE+Qgy7DwfO9iMP2dJgNgltkDyZTrMk+D06QArV2ZvyFnUaIMEZ88G+NnPiASz6L+oOmyNb9k2VXQ8QR9J8KqrANaske3+8pKzQYLXXgtwzz1EgnlaApGgBPq2QOJJEM+H4TkxV5PsdLgcHhDXxdgGCfJeNK0J6vZMuny2xresVqGeIAuisGnTptByRo4cqRVKS1YpUc4WSD6SYDmEg9K1gyxIsJw\/hrq45p3P1viWbVfkxgi\/IyxbmC05GyDxb0ug3uVs\/CwkGLaDSFDNCkVPsJztQK3l7kjbGN8qrYs8IqMbPVqlcllZGyARCcqi746cDU+Q\/4CU+8fQnZ5S08TG+FbRIPGcoEphtmRtgOQrCZbD\/VddOyESpMPSurYTly8yvH4eR2GiFLVBgliXL2uCfDuIBNWGCY8deYJq2JmStjW+ZfULXRPMWylReRv68BfnsT7XHy2P84LEgez6gV9Z4xTlbHiCRIK6vWEun43xraJd5HQYX36LSj7sDvs0HSYSVDH5UlkRu3XrAOrq9MuzndPGh8C2zknlO0eCSQrn8XcbIPlMgq7vcOrakGkCEG0A9XL9nKVpDHT7wmQ+G+NbRb+KvTEiDoBynQ6HDWQiQbkhEIad60eMiATl+lZFKpIEMQI0xhE8cOBAsTzdN0ZUFAqTtfGl8MUTJBJMZ11sOowbSoglJpfXVIkE0\/V3WO7YjZElS5bAxIkTi\/lYhGnVV6HEGyhidGr2EFMUyRIJRh+NIBLUHxSPPgpw6aWF\/DgNZtG5iQT1MdXJaWN8q+gReURm7NixJQTIChXfIZapDPNgQkKNe9cYN2OQEFnsQla2DZDE3WHXp5BRHgCRoIwFhsuw2yLoBe7cCXDyyQU5l22BPEH9\/o7KqXxY2kRkaZ5Ily9fLv3kpqoHGgeXSB5\/\/jPAmWeaB9hUiVHG\/+STAOefX1qLy55MGjxMEwC7LYIkiGvC+DORYJoe0strw8lR0SQTT1BUiL1TMmvWLJg\/fz7U1dUBBmoVvUTRE2xsbIQJEyYYeXf42WcBzj239\/3eH\/3oKNx++wcq2GUqiwSAb7QOGjQIampqinVfeGF1MaAq++W8ee\/DzTf\/J1P9sqgsCgPdumtqCv1\/zjk9sGpVD4wYUfj5xRcPOxtRyDQGutiZyoc2vXv3bmhqanLr3WHGzKbWBHnAeE8Sf4+RqzFIqwwJovzUqVNhypQpqfuguXkA3Hvv\/xTLueyyd6GlpTt1ubYKOHLkCHR3d8PAgQOhf\/\/+xWrq6gbDrl3VUFvbA11dhUfEXW+LLkZRGOiWd\/rpBdfv5ZcL72yyn6+++l+wcOFbusVazWcaA6vKShTe3t4Oa9euDSRNzvQkqi6KZLo7jAQ4Z84caGtrg+HDhxc9P1lPsKWlBUaNGmXEE7zggirYsaNAGpjwaMTmzYdVsMtUFr3kQ4cOBW2vru71YJk3g\/rjFB\/\/4fQOvRnfUhQGuu1k2L33XgEr9AQRP5dtwTQGutiZyoee4MaNG6G1tdU9EjTVSFZO1FoiTo2HDRsWbJqwsP7MM+TzTp482ShIYggl1+\/cRq2HsSMeP\/whwD\/+UR47nLq2ZXJNcMcOXA4paMLWUKdNA\/jlL0t\/p6urrXwmMbClo2q5zq0JqjZARj5uM4XfJMlyd1gMocQPBpk2ZS0TZvz85g4O5L17AT772YJmLu9w6mJnkgD4nWHEChOPp6v4mcRAtx9M53OWBPnHllijxfN9smCwc4C8PF9WHucEmQf1kY8AfPDhfojLu6pJniDqXg6DWNZmwuRMEgDr\/699DQB32FlivycSTNNTanmdJEFx7Y41Cb22FStWFNf01JqqL20DJGbs3\/wmwG9\/6773FEYAf\/gDwJe\/XNBdJMGXX+498qGPvFs5bZCg+OFjduHqHWKTGLjSuzbGt0rbIo\/IsM0KsTD02vAqXXNzc8lRDZVKVWVtgMSMHc+H4dQIk6tff9QtzPjDpnSuD2LVvuflTRFAYyPA8uW9Hw++jpkzAe66K\/xvaXQ3ldcUBqb0MVGOjfGtolfoYenZs2fDwoULgx1cMZk4LK2iIMraAImf9pTrIVnWBv7SP7\/W6fL0XtUGoj4EOuXw94XZeiBfDvu7i\/gRCer0eHweZU9Q59pcWrVtkiAauuvrQFEEEDZY\/\/QngJEj3fVk0tiCKQJguB09CnDccX01Cvu4pNHbZF5TGJjUKW1ZNsa3ik5KkaV1AyioKBQmaxokjBd71lm9RFEOU8gw4w8jQXHHOC32LuU3QQDjxwP85jfxH4kJEwA6Otz8kJjAwKU+tTXTU2mjVmRpvoIsokybJkFxLY2\/Qxo2PVIB1JasaPxh64Gs7nIgdR2cTBCArJfH5O65B2DGDB1t7eQxgYEdzfRLNT2+VTWpyKCqjPTYWtr06QCrVrn55WcdKhp\/1BEPlP\/+9wHuvtvt9qgaKsqnJYCXXgL49KflcHF1bTUtBjq4285DJCiBsGmQGIFcc02BLMrhfF0UCUYt3ru8uC\/R5aEiaQkgaUNErJTJt7YCzJqlq7XZfGkxMKuNmdJMj29VrUI9QXZQmj3AztYCsXDdA9OqivHypkESCaIc1tF442eRY+Ku+vEPCLm4y6ljD2kI4Be\/ALjuOjkvkOnmojeYBgMdzLPIY3p8q+ocSoL8WcCurq4gzD6SH97pZWGw8P9ZJZMg9fT07gjy5OD6Ohpv\/Ozi\/+WXAzz0UHgvPPwwwHe+ozbos+pP3XrSEICqF8h0lF1D1G2Tar40GKjWlZW8yfGto3NiUFVRwXI\/Jxi1ocAHVHDRc2LGf801g+HxxwtRZJL0ZAMY27Z1q455uJVHlwDSeHR8UIW33wY4\/vh8MdHFIF+t42t3ngTR80MlWch7nBpv2LChTwh8myCbBCnuy+7yOhozfhbzTibqzfe+B3D\/\/f54gzoE8POfA2CEHUwYIaa+Xt1SXVpa0MFAvcXZ5sDxfcklD0NHx2VBXNGsU+R0GBVpaGgIgp4OHTo0uCaHsczwZ1S0XKfDcUTH\/nbttQC4huRSQuM\/\/\/xj8PvfF6JKJ3mB4nROhjRdam+YLjoEoDsNFut3hQh1MHC9Xxm2zzyz0x0S5CPIsBfgamtrg1D4+HOWBIgdaMoTbGgAWL062jNKM22ybWho\/HwAVbzzLJPwAaGzzy5I4pSY3ZOWyeuajCoBmCJA8YOi8hEyjaEqBqbrN1ne9u0A48b1lugUCZpsqImyTJGgzKBwdUqc5kC3y+SuYh8qBCDT1yp1o+xrrwEMHdqbS9YbV60nTl4FA5P1mi6Lt8mqqi449dRxRoMmq+hbUYelZQjOlWkP34kXXwzw2GOF33R1HYFTTul9Y0S2s11sl6zuTE6WAPi2mo4M9MILpa8SZk2EshioYpulPN8\/uEyzfv1OMB05XqU9kSTIT4nZ2UD+eUyVStLKmvAE466Z8fr9+98AH\/944TdZG3gYTvyRHnxMad++npI3RlSwteEdqdSfVjaJAPC9HowDiMnmGijaBQbjZQmv1eH1uixSEgZZ6KBbhxjNnfWRifGtqxPmiwygwB5E2rNnD3R2dpZsjLBD1Gkq5vNmEVmaEcCFF\/Z6VVH6u+I18Ye4kQCffvo1GDx4sDYJ\/upXAOyhPpskYcouxHLiCIAfYFm1TfRocJ0W67aZypEEr74aYM2aXlTE\/nGOBMW3f8XQWaaPyGTxxojOGUBm4N\/6FsCmTTbNOrxsngDZ63EHDx5MRYJY04IFAM3NvXW64O3KohtGAP\/6F8AnPhE9wGTL1pUTH+zCEJx4R9lWKicSvOoqgLa2UiQwoC0eW+KTcyQovvgmkqDpw9JZvDanc+qfX4fDdaARI2yZdd9y\/\/IXABbPln01TRr\/7t0AX\/xib7233gpw443ZtU+3JhEDcXqFBL94sW7p+vnE6TGWhA9ePf+8fplROU3agXntCvfw8SzmU0+Vlh7nnTtHgkmeoMnw+qwu2XeHdR5nTrMzyk93ogJwmjYkcVrHplg2jF+czrkaRoxhzDAYN+7U4kPz+Lespr9JfX3FFfiAuPzgTyov7O827EBHDz4PEh8+SbBsWd+SZPrGORLEZqBSYWuCHR0dsGDBAmNb2aLXKRIwg5SB1NjYCBMmTJB+fP2FFz4Co0Z9LCimvf0IXH75MeX+ZmfzMGNr6\/swffp\/lMuQzcDXxTZBeALAh6oHDRpk9G2X4cOr+hCKiw+347T35JN7H51HXBCjP\/6xp7iRJYuzbbnx4z8KW7b0DVl9xRX\/gdWr309VPZKgDTtQVQqJb8WK6tBLBdgvl17aA3fckVwqtmX37t3Q1NRkjFeSay2VkNodZlnYwemwt0dUK0Z5VU8Q80ydOhWmsNX9hErZFTO2qaCjI+Zh5bCBt27dwWAAmkqiZxOm75EjR6C7uxsGDhwI\/furH5GJ0xWnc2ec8ckSEdRhy5Yu+NjH1D8cpnDp6qqCpqaBsGtXX\/K7995D8JnPHDVVlZVy2tr+G2699cTQsocNex+2bu1SrtemHSQp8+qrx8H8+f\/bpz9YPp1x1t7eDmtxWx\/Qi17nzo2RJDC2bt0KZ511FgwYMCBJNPHvKmuCLS0tMGrUKClPkHlVbFMhUZEEAd5LQ1H80j3wgD4R4tTzG9+Q98Lwg3Ho0KGg7dXVpaSQtm0s\/6OPVsHkyVV9ikPj3rOnJ5PgAevXA6xcCaEDDfV44om\/WMXAFJZiOSNGVAfrZWEJ23XqqT3w+OMAVX3hL8mShR1ghfhhxDeZd+6MtzUcX5s3H9beFUdPcOPGjdDa2uoGCfJxA8PC5jPP7W9\/+5uxAAo2dodtHXFpauq77oFG8OCDve\/\/xg0iPKKCl\/jDFo1\/\/eveB5LCyshqLQgHKp5NDHlosKgWtvkHPwDAnb40accOgJtu6osHXybWhcea8C53VhikaZNMXnFHOSoPth37Ac8gnn56Qco0BniN9IEHChsaUSQt9ge+1Y2R9EwdB3JmTZBfB8TprrgBwhRFQJYsWQITJ06U6W8pGZPnBHkCPHAAYPBgKRWUhL7+dTOhqdCIcOe5phATITaZNv6k+tjf8X1eHbJjA0RmYIm6YN4bbgA4\/\/zSgZYXBrJY6ci9+WYh7qP4YdQpiyclHdzD+gF\/N3cuAAYVsZWcIEHm4fHBEfhNi1dffTXYEDG9JigLqgxI7e24Xthb4htvAAwaJFuDuhwaGRowf9REppSoAZ6U1xUCkPViktrD\/o544L+WlmQsXcFAtm26cn\/\/OwAGzJX1znTr4fsA\/48fHXxuImlKnrY+Mb\/M+DZdJ19esDEihtNHAUaM\/\/znP2H79u0wY8aMzKPHMEWTQHLpqMejjwJs2VKYUuIhXtzDOfPM9F3oOgG8\/nohQs+2bX2vG55xBsBppwHg4dlTTtHHwnUM9Fsmn5NhcNxxg2Hr1mp49lmAvXsBPvggvAw834ofmYsuKtghP1bka7UrmTS+7db+4bU58agKq5RNU01Pf1UbJYKEX0hcPOdvPmCZMmeSVOt2RZ4IwPx6mCt9q6KHj3bgPAkeOHAguDdcI7NwpdKbCrII0tlnR0ec9Zn8GEw+Gr+CCQSihIGfGDhPgmh8WQdRDVszCCNBJL8sLq2rDlYb8kQAfhKAqq34aAdEghJWgCDhgvykSWOCU+gnnSSRyTMRH41ftYsIAz8\/BESCEiMhb5AkVLQuQgTgJwGoGo6PdpD3+C7ZHX7uueek+iTsILVURk2hvEHSVNtoNh+NXxUgwsDPD0He47uiwutl6+AQAAAIeklEQVSrDjqX5IkA\/CQAVRvz0Q6IBCWsIG+QJFS0LuKj8auCRhj4+SHIe3yTJ6g6EnOSJwLwkwBUzclHOyASlLCCvEGSUNG6iI\/GrwoaYeDnhyDv8U2eoOpIzEmeCMBPAlA1Jx\/tgEhQwgryBklCResiPhq\/KmiEgZ8fgrzHN3mCqiMxJ3kiAD8JQNWcfLQDIkEJK8gbJAkVrYv4aPyqoBEGfn4I8h7f5AmqjsSc5IkA\/CQAVXPy0Q6IBCWsIG+QJFS0LuKj8auCRhj4+SHIe3xn4gmyAK2bNm0K7H78+PEl4blMhtdXHVjlIk8E4CcBqNqfj3ZQESSIjylhwndJ4h53x7vLSIhr1qwpeckub5BUDdWGvI\/Gr4oTYeDnhyDv8Z2JJygaO\/\/C3PLly2HYsGEBQUZFuM4bJNXBakN+\/\/79wfus+O4y4lWJiTDAd0f8s4O8x3cuJIjeHqZZs2bB\/Pnzoa6uLnh0WfQS2UBnIDU2NsLo0aMrcfzD66+\/Dk1NTUAYEAa+2QGz7bJ6fD0NCyGhsSkvljNt2rQgcnUcCXZ1dQUEsGvXrjRVU15CgBBwFAF0blpaWqC2tjZzDY17guImCP9IExLgnDlzoK2tDfBtYyab5AkiKkiE+I8SIUAI+IcAkl8eBIhIGifBqO7hPcABAwYUxdArTFoT9K\/LqUWEACHgCgKZkGAUASII\/CZJ1O6wK2CRHoQAIeAfApmQIDsHyMPHnxVMOifoH+zUIkKAEHAFgUxI0JXGkh6EACFACIgIEAmSTRAChEBFI+A8Ce7btw\/q6+vhwIEDMGPGjNwfgrdlLbg2umDBgmLxfFvjMPBlKYHfIGMgRLWNHarHNWTxCma524uIA98exIV\/6dFnHGyNs7BynSZB\/gYJdj4erB47dmxwu8S3FEYC2MY4DHzZVGJkxx+nimpbTU1NiR2wg\/d41rTc7SUMB9xUXL9+fclde7QL8WKBTzhkPbadJkH8Ci5evBh++tOfBneJ+YGBg8GXJJ6X5NsVh4HMlUOXMWLtHjJkSKAmOyqF\/486OvWpT30KZs+eDQsXLgzOmvInD958882ytJc4HNDm8aockjyfkPB9wyEvW3WaBMWvYNxRm7wANFEvP60RpzxRGNx5552wdOnSxCuHJvTLogye9OIO0X\/hC18oITr8SMydOzfAAkmQ95rK0V7EGYF4soJdLRM\/jr7hkIXNsTqcJkHR8ytHo5bpTN6A0btBw8c10ObmZujo6IDOzs7idIhhgFeM8Cph0pVDmfpdkIk7NM9P\/U477bSSSEM8dnv27AnFSoxK5EJ7o3QI+xiwJSD+xhUSPh9xyTccsuwjp0mwUjxBscNlvurkCRaWSGSwKlcSFO2C\/xj47hETCX6IQKWsCYaRIFsLjVvnKvc1Qb7dYdPAsOuUPq4JxuHA\/41fJvAdByLBDxEo990+2Y4UPV7ZnT5fdofFjRD8uRJ3h8M2hKI2P3zeJZcdN6bknJ4OYyPL\/dyXbEfx5wRVzr7ROcHSpxrK3V7izgniLjqLwIR2RecEZUdXvJzzJGimmVQKIUAIEALhCBAJkmUQAoRARSNAJFjR3U+NJwQIASJBsgFCgBCoaASIBCu6+6nxhAAhQCTogQ2EBa1lzeKjjnjQVKtN4O\/piu\/h8BXzd33FO71hCor3fK02ggpXRoBIUBky9zLw1+x8CiyRJdLiwXyTJIjtiIoGk2UbqS7aHfbWBogE03Vt2HvXpkkwLlJQOu0pd1oEyBNMi6AD+WVIkE3Jxo0bB7feeivwB2\/FgK5hj2DzU24M+IpX2lhgBxzg\/PvRCEnUlJEvRzz8y4JDNDQ0wMyZM4vIivpERd3BDKjHpEmTSmJOxhEa5hEDWDDPjX8eNm46LOLHZMX2+RoKzoEhkEoFIsFU8LmRWZYEkSAw8QEFcGCuWLGieBOB3bi4\/vrri0SC5SORsHyMyNjNFlkSFMvBMidPngyM5NjP4iNcfN1h+mEbNmzYEOi3atWqYgQetjSQhE8YOaXxBBlJjxkzpiQOYBjZumFBla0FkaAH\/Z80yLGJbGDyXlLY71CWJxUM4IDPGyxbtgxwUPNeHv4fw33JkCAjL74czM\/rjuHyRe9LzJfkTYnkFdVG1u1RHisj5DjzCHvugZXHsOHXaPm78AxLD8yv7JtAJFj2XVggkpUrV\/ZpCT8dCxuAUd6OGJ6Kj1vHKuHJS4YEeWLFEFgs8TEiX3rppZIYeWyqykiYPbGA7YralRVJL2kqHEWSup6g6FmHTaN9fSKiXIcSkWC59hynt4onyIKwYvY4b4cRKAYqZVNNnrx4j0yWBPmHpHjY2TEeWRJMIhFeNww3xgLUhu2cmyRBcXovmpbq0RoPTLMsmkAkWBbdFK9kGhKMWvwP89R4EjTlCfItC4sczk+HZTxB5j1iyH0kfNRT3Cjh6zRFglHrgOQJuj\/AiATd76NEDXVJMGqdTny8SHZNMGy9kW0ORE0vk2IiijrKtJV5XMcffzzs3bs3eH8Eny0IS3FrgrK7w1iujF60JphoyrkIEAnmArvZStMMQHENK8wzEtfz2JEQtouLrcHnUDHhRglOO9k6Jb95IO4OiwSX5AkioYbtDocRrKhj3CHytLvDceuAfE\/T7rBZuzdVGpGgKSRzLCcNCaLa4jk3\/v1f1iwx6Ct6We+8806R9MSze1gGPhWJid\/EEDdx+DOAMiSI5Yl1iefx2JQYPVj+qE9UF6U5J4hnGvHoEe5shyW+fUk72zmaUEVXTSRY0d2v33gZ4tUvPX3OpF1hvoawGyPpNSgtgW6MmEbUXHlEguawrKiSXCdB1E\/0QuM6SLw7bLoz6e6waUTNlUckaA7LiirJVRJUWQsUO4yPImOyMymKjEk0zZdFJGgeUyqRECAEygiB\/wNXxAIiISaukwAAAABJRU5ErkJggg==","height":0,"width":0}}
%---
%[output:22b1053c]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAUEAAADCCAYAAADXXXDzAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnX1wFlf1xw8tPxK0yBjNRCi1UE1nin+gZrShKFRq25GaIkwViJ2EkAawBWILKRCFWqmSEugYykt5iTGZkVCxTqHWmVo77YiTwjgojNJWMi2Z\/pCA+RVtC1NCafnN2Xgf7rPZl3ufZ3fv7tmzMwyE3L177\/ee+9lzX\/aeIZcuXboEfLECrAArkFIFhjAEU9ryXG1WgBWwFGAIsiGwAqxAqhVgCKa6+dUqf+DAAaisrIRdu3ZBeXm52k0JS9Xd3Q01NTWwePFimDVrVsJKz8XNRwGGYD7qpeRehmBKGjql1WQIprThdarNENRRi9MmTQGGYNJazEB5GYIGROdHRqYAQzAyqZP7ICcIHjx4EJYuXQpf\/vKX4Uc\/+hF87GMfc6zge++9B52dnfDLX\/4Sjh8\/Dp\/4xCdg3rx5UFVVBR\/96Ecz9+BOraNHj8KWLVvghRdeAPy5rKwM7r33Xpg0aRJcccUVVtpz587B7t27rflJkd\/s2bOt\/IqLi600Z86cgdraWrjrrrusf2\/atAmuueYaePzxx+GGG26Avr4+2Lp1Kzz11FNw\/vx5uO2222DGjBmwevVqnhNMrpnmXHKGYM7SpedGOwT\/+te\/wrJly2DChAmeADx79iw89NBD8Mwzz8Cdd94J06ZNg1deeQXa29th8uTJ8PDDD8NVV11lCfmnP\/0J7rvvPrjuuuuguroaCgsLYc+ePbB\/\/35oamqCmTNnwr\/+9S\/ruQhgXMSYOHEi\/O1vf7PyGzNmDDz22GMwbty4DATfeecd62eE5MmTJ6GiogLef\/99WLJkCbz++uvWc8aPHw+\/+93v4LnnnoP+\/n5Ys2YNL4ykx7StmjIEU9bguVRXhuDIkSNh0aJF8PnPfz4LYk75Ivzq6+th7dq18J3vfAeGDBmSBTz0JNGD+\/e\/\/20BEIG4YcMGGDFihJUOIbpq1Sq48sorLS\/t6aefhubmZsuzmzJlSuaR\/\/jHPyzP7xvf+AYsX74cEH74M\/7d2toKY8eOzaTt6OiwnrF582b4yle+Yv3\/xYsXLS8R\/2BZeXU4FytJ7j0MweS2XWQlFxD84Q9\/aA1r0XtCWAgvzqkg6FU1NjbCq6++aoFo1KhRmWQ4RF6xYgXg3+i94fYUhCF6fHfccYdjvRCIDzzwgPU7vEd+NkIMPbg\/\/\/nP1rMKCgosCH7605+28hw+fHhmKI3PxQv\/Xx6OHzt2zPIu0UtkCEZmWrF4EEMwFs0Q70IICIpS4hBz+\/bt8JnPfMa14OjdzZ8\/H4YOHQp1dXUwbNiwrLRPPvkkvPnmm7Bz505raOq3D1HM8+E+RfT27Bfmh55cW1ubNe+IEPzSl74EK1euzHigXnmI3+HQmSEYb3sMunQMwaAVJZifgCAOaREQOC\/3hS98wfK+hJdlr7aAypEjR1wVGT16tAWtt956K28I4mIJDpNlCNqByRAkaJwBVIkhGICI1LOQ5wRvvPFG+NWvfmXN1YkFC6w\/ptm2bZslBQ5v8d8\/\/\/nP4cKFC4OGnna9cKHFbTiMefzhD3+wFmDWr1+vNRy2Q9A+DJeH1PzFCHUrdq8fQzC9ba9cc\/vq8Lvvvmttj+np6bGGxbjw4LSNxmkRAh+Kq7w494ZbWn7yk5\/ABx98YC2M4DBWnms8deqUtWUFh9\/odeKQ12thZOrUqYDzljh\/iMNhp6EzlsmeB27Hwf\/H1WpeGFE2CzIJGYJkmjK8ijgBDr23uXPnwje\/+U1r5RZ\/fvTRR62FiaKiIqswYovMb3\/728yWlt7eXgtmb7zxRtYKrX2LDN6PW19OnDhhgRaH315bZD7+8Y9bw+Hrr78+s0XGCYKiTOhd4lwlb5EJz26SkjNDMCktZbCcThDEFVncZoJ\/fvazn1ngs0MQiyw2N+PG5Ndee83a\/vLVr34VFi5cCJ\/73Ocyixb2zdJ47y233GJ5nPICjM5mabdFFMwDPT8car\/99tuWx4grwxs3brT2FPLCiEFjM\/BohqAB0Sk+EkHpBEGKdeU60VKAIUirPY3VhiFoTHp+cJ4KMATzFJBvH1CAIciWkFQFGIJJbTkuNyvACgSiAEMwEBk5E1aAFUiqAgzBpLYcl5sVYAUCUYAhGIiMnAkrwAokVYFIICg+V8KjlfDCc93E6R7icyU87w0vPKNO3nCbVGG53KwAK5AMBSKBIH4hgBduQhVAvOmmm6yfcVURTx6WjzxKhnRcSlaAFaCgQCQQtAuFUOzq6rLAt2\/fPusbVKfjkSgIzHVgBViBeCtgBIL4ZQFeCD78tzh9BP+PcmzbeJsCl44VSKcCkUNQ3lSLZ9HhSb\/y0BjPqsMz4UpLS7NaBD+kxz98sQKsAD0FMEYM\/jFxRQpBBKAb5LDy9vlCIQjCr6GhwQqwwxcrwArQUwDPqcQjzkyAMDIIqnxWJSA4Z84c62QPcYlTTFCkq6++mp4FKNQIXwAtLS2WobAGrAElOxC2bWoqLBIIugEQjzvHo5IwIA8Of93SpSH4tx8HcfEIz9fDMJFy9DS\/+yj9njUAaxGRmh2Y7t+RQNC++IEdU+wVxKEunuWG+wRFzAn7fKBpkeIAEgwSjgeSYtQ2jMmbxos1ACtYPDU7MN2\/I4Fgvh3WtEj5lj+I+ykav64urAFDUNdmVNIzBFVUikEaBgBNAOiaFkU7MO3kMAR1rdBQeorGrysla0DzRcAQVOgJpkVSKGLoSRgANAGgazgU7cB0\/2ZPUNcKDaWnaPy6UrIGNF8EDEGFnmBaJIUihp6EAUATALqGQ9EOTPdv9gR1rdBQeorGrysla0DzRcAQVOgJpkVSKGLoSRgANAGgazgU7cB0\/2ZPUNcKDaWnaPy6UrIGNF8EDEGFnmBaJIUihp6EAUATALqGQ9EOTPdv9gR1rdBQeorGrysla0DzRcAQVOgJpkVSKGLoSRgANAGgazgU7QD796xZi+DJJzdlnR6lq02u6dkTzFW5iO+jaPy6ErIGdF4E998P8PTTeCrOZSt4+eUDDEG3TsGeIB3j1wWfnJ4hmEw7uHQJYOpUgJde8m59hqCHPgzBZBp\/PsBzupchmAw72L8fYPVqf+iNHQtQXw9QXn4AKisrjcUX4uFw0D01pPwYAMkAQEjNn8k2jnawZAkAhhSXh7Z2HRB4+Oe55wCGDcv+rWknhyEYttUGlH8cjT+gqilnwxrE40Xwta+peXlf\/zrAjh3+zZsKCIrYIc\/g60I6VRqjzeElTp7mk6XdDYYBEA8A+HfpcFNEbQcXLgDcfrsa9FatApg3T7\/+qYAgBlvHa9asWYMiysmB2I8cOWIBsbW1FYqKijJqmhZJv1mDvyNq4w++BvnnyBqE\/yJ45RWA++7zhh4Oa\/E6fBhg5Mj829V0\/zYyHJbBt3HjRitwEAISAy\/V1tZaQdmdos2ZikaVfzPnnwMDIHwA5N9K4ecQtB3s3g2wbZs\/9BB8L74YTv1SCUH09vBasmSJFXxdhNh0izssRKqvr4fp06fDpz71qXBaI8a5ovGfOnUKSkpKQEwjxLi4oRSNNRh4EeRqB7hw8ZvfDIPnn7\/Cc7vKmDEXYcqUi7BzZyjNmJUp1uXQoUNWXHFTTk7knqAcVhPVkD0\/Pwhiegw5WVVVFX7rxOwJ\/f390NfXB8XFxVBQUBCz0kVTHNYAQFeDpqYi2L7de8yK0Js27Rzcffc7gP+O8uro6LBCiOKVCggiAJctWwZtbW1WnGF7sHU\/CGLA6bKyslR6gqjN6dOnrbqnNeQmawBWn\/GygwceANi61TskKw5t77\/\/Q5g\/\/0KUvHN8FnqCe\/fuhZaWFvoQlD1AedEDh8Y8J+hvi0HPBfk\/MX4pWIPsedFTpwrhiSdwd4V3WyH0NmwAmDkzfm2KJUrFnKAbAFEAXh1WM0wGAC+MoKWsWfMhPPvsBTh40N3bQ+ghHHFrSxKuVEBQ7AOUG6SiogKampqsSX7eJ+hvqgzBdEJw166BDcde390mDXp2a08FBP27uHcK0yLlW\/4g7mcI0ocgrt7i5uQFC7yhh4sXTU2XYNKk\/7E+RUv6Zbp\/R746nEuDmRYplzIHfQ9DkCYEP\/gAAD8v8\/P0Fi4EWL6cpgam+zdDMGhahZQfQ5AOAPy+vUXv7s47AVpaBhsTRTtgCCpAw7RICkUMPQlF49cVLaka4GdoW7a41xahd\/PNAG1t\/ookVQOvmpnu3+wJ+ttdLFJQNH5dYZOgAc7rvfEGruK6D3HFsVIIPd05vSRooNuuDEEFxUyLpFDE0JNQNH5d0eKsgcoQ98ABgJIS3Vpnp4+zBrnWzHT\/Zk8w15aL+D6Kxq8rYVw0QG+vudl\/iHvPPQDf\/a6+t+elS1w00G07Hg7nqZjpN0WexQ\/kdorGryuMSQ3efhvgW98KZ4iro4NJDXTKqZPWdP9mT1CntQympWj8unJGqQF6ezU1\/ltXnn8e4LOf1a1J7umj1CD3UurdyRBU0Mu0SApFDD0JRePXFS1sDd58E08p8vb28Ptb\/A7X1BW2BibqZbp\/sydootVzeCZF49eVIWgN0Nv78Y+9t6aIw0R1V3F166aaPmgNVJ8bZjqGoIK6pkVSKGLoSSgav65oQWjgN8xF2NXWAtx9d7ALGrp1dUsfhAZBlSWofEz3b\/YEg2rJkPOhaPy6kuWiAUJv61aAdevcnxY3b89Ll1w00NU56vQMQQXFTYukUMTQk1A0fl3RVDVQ8fZmzMDwDvH09hiCupaRX3r2BPPTL7K7VQEQWYEMPMhLg7NnASoqvBc1XngB4LrrDBQ8wEdStAPTTg5DMEADDTMrisavq5esAZ6qvGcPwIMPOueSz6dpuuWKMj1FO2AIKliQaZEUihh6EorGryvaa6+dh9Wr34c9e0Y43orgq6sDqKxM3jBXVQuKdmC6f0fuCcoxRbDhu7u7oaamBk6ePGnZwYQJEzj4ukOPoGj8qh3f67tcBN\/vfw9QWqqaW7LTUbSDVEFQHKO\/du1aK9g6XihAZ2dn5qh9JxM1LVIcug1F43fTVWVhAwOBx2XvXpT2QdEOTPfvSDxBEUpz9OjRlr2I6HL4bwy01NPTA8vx2FyXy7RIURq527MoGr9cVy\/wifm9rVvPQ0FBL4waNSq1YUcp2oHp\/h0JBGVjtw+H7UGYnAIwC5Hq6+th+vTpqYw7jMaPMVpLSkqs4FQUru5uPEG5EBCA9gvjaNx220VoaLjs8VHUQLcdqWmANn3o0CFoaGigH3dYNLYMQXuwdYSdHJxd3CMgiD9XV1dDVVWVru0kPn1\/fz\/09fVBcXExFBQUJLY+x44Ng3vuKYETJ4Y6gq+6+h24\/fZzgBC0X1Q0yKfxqGnQ0dEB7e3tliRODlA+Wqnea9wTlAtqh6Idgs3NzVBWVpZKTxC1OX36tFX3wkL3mLOqDR9lur\/\/HWDGjKGO4MOh7vr1F+GOOwZDz17GJGsQlN7UNEBPcO\/evdDS0sIQRCMREJwzZw6Ul5dn7Mb0nEFQBpxPPkmbC3r99YEoak5DXQTfL34BMGWKniJJ00CvdmqpKWpgun8b9QTPnDkDS5cuhcbGRigtLbVWinG43NraCkVFRQxBqV8kwfj9Fjdw1DN5slpnd0qVBA1yr53anRQ1SDUEsdnlfYK4etzW1mYBUb5Mi6RmnuGmiqvx+4HviScAbr89GG3iqkEwtVPLhaIGpvt35J6gWlNnpzItUi5lDvqeOBm\/H\/gaGwe+3Aj6ipMGQddNNT+KGpju3wxBVesznM608SP4HnkEoLV1sBA4xxfFicumNTBsAtbjKWrAEFSwLNMiKRQx9CSmjP\/Xvwb49redwYcBwx96KLovN0xpEHrjajyAogam+zd7ghoGaDJplMaPXh9+r+u2snv8uBklotTATA39n0pRg9hDEFdwDx8+DFOnTnVsIdzWsnHjRqirq8ta0fVvTvUUpkVSL2l4KaMwfreDCuJy8nIUGoTXgsHkTFED0\/3b1xNECNbW1sLs2bMzhx6I5hSFdzr5JZgmH8jFtEhB1iXXvMIy\/m3bABYudB7ubt4MMG1ariUO\/r6wNAi+pOHlSFED0\/3bF4IyhMTpLwKMR44cgQULFngefhCEOZgWKYg65JtHkMbvtrob94NIg9Qg3\/YwdT9FDUz3byUIyiAUjR+29ycbmWmRTBm8\/NwgjH\/+fIAdO5y9viQcTRWEBnFoy3zKQFED0\/1bGYJOHmE+jalzr2mRdMoaVtpcjd9tkQO9vlWrAObNC6vEweebqwbBl8RcjhQ1MN2\/tSBoan7OtEjmTP7yk3WN\/5573Pf0mVrdzVdHXQ3yfV4c76eogen+7QhBec5PxRDCHhqbFklFg7DTqBr\/uHGDt7ag1\/fsswDjx4ddynDzV9Ug3FKYzZ2iBqb7t7YnaMIETItkos72Z3oZf0sLwPe\/7zzXl1Svz0lzigDQtS2KGpju3wxBXSs0lN7J+J329aHXt307wK23GipoiI+lCABduShqwBBUsALTIikUMfQksvHfcMPgI+kRfpS8PvYEnU2KIRh8V2NPMHhNQ8mxvf0DmDv3yqy8EXx4cKnTtpdQCmE4U4oA0JWUogamnRyGoK4VRpzebcibpli7QnKKANA1J4oapA6C9mhzaAQi4hwfqnq5SzjBD4MPdXdfTFyMEd2O7paeIgB0taGoQaog6BR8HeMOd3V1WcHX8TO8tB+v77bF5dVXz0NvL8fcZQ3o2UEqIOgVfF32DMX+RAzEnrZAS0OGZPsEON83d+7AeX14UfQA2AvSVYCmHaQCgnJTO8UdFtHl\/EJumopLqm+qane89x7ARz4yGH779wOMGZP9\/wxBmgBQs5TLqSjaQaohaPf8\/CBYX18P06dPT3zc4VOnhsC4cdkB1MV8n9d8GMZoLSkpgeHDh+v2HRLpEQCsAS0NsD0PHToEDQ0N6Yw7bI8z7AdB7MnV1dVQVVWVyE7d2zsUJk26JqvsCL8\/\/vF\/fevT398PfX19UFxcDAUF2QD1vZlIAtYAgJoGHR0d0I6xWAHSCUGxMjx27FjrwFa\/OcHm5mYoKytLnCf4n\/8AjBpVmIUinPPDxQ7VC18Qp0+ftupeWJidl2oeSU\/HGgBQ0wA9wb1790JLS0t6IUh5dfjSJYArrhg855fLlx0U54J0ocwa0JwXTfWcoOgEFPcJOq325gI\/oREDgCYA+EVgPnwGfzGia4U+6e37\/IL6ppchyBBE06NoB6nzBHNhjmmRVMps\/8IjKPixJ0h7e4iKbclpGIK6ivmnZ0\/QXyPPFHV1ADt3Xk6C8HvpJYBrr80zY9vtFI1fVyHWgD1BXZtRSc8QVFHJIc25cwBXXZX9i\/vvB3jssRwz9LmNAUATALrWQtEOTI\/0GIK6VggAYc37eRWFovHrSs8a0HwRMAQVeoJpkUQRw573Ywh6GwNDkCGogAvtJOwJKkj24YcAV0rnmeK83+HDACNHKtwcUBIGAE0A6JoHRTsw7eQwBH2s0MTQ16lIFI2fAaCrAM0XAUNQwQ5MiHTwIEB5efaqbz6bnRWq6ZmEIUgTALp2QdEOTPRvWXf2BB2s0O79PfUUwMyZuuYabHqKxq+rEGtA80XAEFToCVGJdPIkwNVXx8f7k6VhANAEgIL5ZyWhaAdR9W83rdkT\/K8ydu+vqwtg4kRdEw0vPUXj11WLNaD5ImAIKvSEsEWSDzsI+nM3heopJWEA0ASAUuNLiSjaQdj920\/jVHuCjz4KsGLFZYkqKgD27fOTzMzvKRq\/rpKsAc0XAUNQoSeEIZI8\/I2r98dzgtnGwRBkCCrgQjtJKj3BJAx\/7S3JAKAJAN0eS9EOwnBydHRNFQTtJz3Png3Q2akjl7m0FI1fV03WgOaLIPUQ7O7uhpqaGjiJ+1MAYMKECdDa2gpFRUWZPhKESOvXAzQ0DGSZhOEve4KDEckQZAjqvjhV0hv3BBFwnZ2d0NTU5BpKMl8IygcfJBGA2JAMANaAqh3k279VQOeVxjgEMdBST08PLF++3LWc+YiUtAUQNxEYggxBhmC+uHO+3zgERZAlUbxdu3ZBufzRLuQeiIUKAKkav65J84uA5osgHydH14ac0huFoD3YOoqxbNkyaGtrg9LS0kFzgvX19TB9+nTfuMMnTgyFW28dCj09A1lggPPu7otB6GUsDwQAxmgtKSlxnTYwVriIHswaDECQkh1gXQ4dOgQNDQ3piTvs1V\/sUBRpxZsCf66uroaqqirXbBCAlZWjAP\/Ga8qU96Ct7VRE3TS8x\/T390NfXx8UFxdDQUFBeA+Kcc6sAQA1DTo6OqC9vd2yOqdRYBTmaNQTtFdQQHDOnDlZQ2IBwebmZigrK\/P0BEtLh2YAuHjxJVi3rj8KHUN\/Bmpz+vRpq+6FhYWhPy+OD2ANAKhpgJ7g3r17oaWlJZ0QPHPmDCxduhQaGxut4S\/CDucIc90iI88B3nsvwObNcezKuZWJ58NozofpWgNFO0j1nCAagLxPcPTo0YPmAzGNikgUtsF4dQiKxs8A0FWA5otApX\/rK6V+R6yGw27F9hPpxRcBpk4duDup+wD9mowhSBMAfu1u\/z1FO\/Dr37oa6aYnAUHxLTBVAGKjUjR+XWNlDWjaAUNQoSd4iZQGADIEB4yEIUhTA4ZgHhCU5wHxcATKFwOAJgB0bZaiHTAEFazASaTeXoDRowdu\/sEPAB55RCGjBCehaPy6zcEa0HwRMAQVeoKTSGkZBgt5GAA0AaBg\/llJKNoBQ1DBCuwipWkYzBC8bCAUAaBg\/gxBXZE00ydydVh4gdOmATz7rGaNE5qcAcCeINXFIfYEFaAki7RyZTm89BLd\/YBucjAEGYIMQQVY5JAkcZ7gxInlVjWPHgUYPz6HGif0FoYgQ5AhGE7nTRQEe3t3wfnz5WS\/CvFqYoYgQ5AhyBCE48ffsFTYsgXge98LR5C45soQZAgyBMPpnYnxBG+99QScPXuXpQL1jdFOTc0QZAgyBFMOQTEXSPn7YB4Oexs5vwhovgh4dVgB7iiSgGAavUCqHoBC02clYQgyBHVtRiV9YobDDMHz0NvbC6NGjUrtydIMQYagCtR00yQKgmkdCrMnOGDWDEGaGvBwGMA6Un\/btm3gdbI0eoI33wyAB6im8WIA0ASAri1TtIPUQxCDr3d1dUFTUxMcOXLENcYIQrC9\/f+gquqTunZDIj1F49dtGNaA5osg9RBEL3Ds2LEwa9YswMBLtbW1sHz58kHR5hCCL798YFBgdt2OlNT0PT09VmhCDDmKeqXxYg0AKGqQagjaQ2x6xR1GCK5f\/xRMmnR1Gvs\/\/POf\/7QCVGMA+htvvJE1YA3I2IGw7VTGHbZ7fm4QPHHiBODxWefONUBh4cFUAoArzQpQVgBf7BhXfMyYMZFX0+jqsKoniKogCPEPX6wAK0BPAYSfCQCikkYhiAVQmROk1+RcI1aAFYiLAsYhqLI6HBexuBysACtATwHjEBTeoNc+QXqyc41YAVYgLgrEAoJxEYPLwQqwAulTgCGYvjbnGrMCrICkQOwh2N3dDTU1NXDy5ElYsGCBtZGa4oVzoytXrsxUTa6rlwZ+nxwmRSt5gUyU2a1uYmsVfmFUUVFhfW00fPhw67ak24tdB7k+WL8JEyZAa2srFBUVZT4uoKhDlHYbawjK+wix8VesWAE33XST9XUJtcsJAlhHLw2oLCoJ2K1duzbTtm51Q9jJdoD34oUvx6Tbi5MO+DVFZ2dnFuixvvY9tZR0iLpvxxqC+Bb86U9\/Chs2bLDefHLHEG\/+qAUL43n2\/ZLyM7w02Lhxo+8nh2GUN6g8Rb3x4Ay8xOeT+G+3rVPXX389LF26FBobG6G0tBQQEpgWvaO33norkfbipQPaPH4qZx8BIfCp6RCUXenmE2sI2t+CssEjFKlc8vDOPuRx02DTpk2wbt06mDNnjvU9tdvXNknRSIae1yb6L37xi1mgw5fEgw8+aGmBEJS9piTai31EILxD0Y7i0zL7y5GaDlHabawhaPf8kmjUKo0pGzB6N2j4OAeKc1379u3LnLKD3q\/QAD8xwm+JxWETlCDo9Tnltddem3XSkKzdX\/7yF0etxByaSluYTuP0MhBTQNj2y5Ytg7a2Ngv4wgNGh4CaDlG2Q6whmBZP0N7gKm919gQHpkhUtEoqBO12Ib\/oqHvEDMH\/KpCWOUEnCIq5UK95rqTPCcr1dhoGOh2xRnFO0EsH+XfyNAF1HRiC\/1Ug6at9qg1p93hVV\/qorA7bF0Lw5zSuDjstCLktflBeJVftN0Gli\/VwGCuZ9H1fqg0l7xPU2fvG+wTTs0\/QHn6C8n5J1X4TRLrYQzCISnIerAArwAq4KcAQZNtgBViBVCvAEEx183PlWQFWgCHINsAKsAKpVoAhmOrm58qzAqwAQ5CADdg\/rZKrJJ86QqCqoVZB\/k5X\/joDv+KRL\/lbX5VTjezf+YZaCc5cWwGGoLZk8btB\/syO0sESUSpt35gfJASxHm6nwURZR36WswIMQQKWwRDMrxGdvrsOGoJeJwXlV3q+O18FGIL5KhiD+1UgKIZkkydPhjVr1oC88dZ+oKtTEGx5yI0HvuInbV1dXdYhD9jBa2trM4c5oCRuQ0Y5H\/vmX3E4RF1dHSxatCijrL08bqfu4A1YjtmzZ2edOekFNLzHfoCF8NzEYQV+w2G7fqLg9vpRPQouBl0gryIwBPOSLx43q0IQAYGXfKAAdszHH3\/cOpkEO7v4Qmfx4sUZkGD+CBJxnwCZ+LJFFYL2fDDPyspKEJATP8tfzNjvcSof1mH37t1W+Xbs2JE5gUdMDfjp4wSnfDxBAWk84kyeM3SCbTwsKN2lYAgSaH+\/To5VFB1T9pKc\/g\/TylDBAxwwvMH69eutcwtlLw\/\/reoJCnjJ+eD9ctnxmHi792W\/z8+bssPLrY6i2d08VgFkL\/NwCvcg8hPayHO09iPCCJgeiSowBAk0o9vqsDwcc+qAbt6O\/Xgq+dw6IZcMLxVPUAanZHMDAAACCUlEQVSrfCCufEbksWPHss7IE0NVAWERYgHr5bYqa4ee31DYDZK5eoJ2z1o2r6Sf+UigqzhWgSFIoGV1PEFxCKuY98LhqNMlAIoHlYqhpgwv2SNThaAcSEp+ptjGowpBvzgzctnwuDFxQK3TynmQELQP7+266m6tIWCaiagCQzARzeRdyHwg6Db5L57odpp3UJ6gXDOnZ8nDYRVPUHiPeOQ+Ah\/LaV8okZ8ZFATd5gHZE4x\/B2MIxr+NfEuYKwTd5unswYtU5wSd5hvF4oDb8NLvTER7GVXqKjyuESNGwNGjR634I\/YVXpU5QdXVYfvcptteTZ4T9DVlIwkYgkZkD\/ahKmBw64D2OSwnz8g+nye2hIhVXKwNhsHES8QAFvOU8uKB20qvWCzx8wQRqE6rw06AtZfRaxN5vqvDXvOAckvz6nCwdh9UbgzBoJQ0mE8+EMRi2\/e5yfF\/RbXsh76il\/Xuu+9moGffu4d5YKhIvORFDPsijrwHUAWCmJ\/9Wfb9eGJIjB6svNXHrYny2SeIexpx6xGubDtdcv38VrYNmlCqH80QTHXz5155FfDmnnv+d\/qtCstPiGLVlr8Yyb9Nw8qBIRiWssTzjTsE5TgtKk1h\/3ZY5R6dNPztsI5a0aZlCEarN5mnxRWCOnOB9saQT5EJsqH4FJkg1Qw+L4Zg8JpyjqwAK5AgBf4f9fTOqn4C9RoAAAAASUVORK5CYII=","height":0,"width":0}}
%---
