% theta_correction.m
% 功能：theta角度校正系统的完整实现

function theta_correction()
    % 加载实验数据
    global DATA;  
    DATA = [
        148, -7.35, 0.12;
        139, -0.56, 5.86;
        80, 0.27, 10.62;
        235, 0.5, 15.27;
        287, -0.35, 20.13;
        253, 7.43, 25.23;
        245, 12.63, 30.65;
        287, 10.84, 35.45;
        258, 17.57, 40.06;
        308, 10.84, 45.02;
        205, 33.46, 54.8;
        208, 33.93, 58.88;
        224, 23.34, 62.6;
        219, 23.66, 63.36;
        215, 17.99, 65.75;
        214, 15.74, 66.32;
        209, 13.01, 67.4;
        204, 10.42, 68.38;
        196, 8.48, 69.55;
        184, 9.44, 70.83;
        183, 1.44, 71.53;
        166, 13.75, 72.49;
        161, 12.91, 73.15;
        155, 0.8, 74.45;
        144, 3.8, 75.46;
        136, -4.91, 76.5;
        126, -7.44, 77.57;
        115, -8.3, 78.59;
        104, -13.32, 79.73;
        96, -15.1, 80.5
    ];
    
    % 提取数据
    r = DATA(:,1);        % 距离（米）
    alpha = DATA(:,2);    % 仰角（度）
    theta = DATA(:,3);    % 方位角（度）
    
    %% 第一步：多项式拟合
    fprintf('\n=== 步骤1：多项式拟合 ===\n');
    
    % 将角度转换为弧度
    alpha_rad = deg2rad(alpha);
    theta_rad = deg2rad(theta);
    
    % 转换为笛卡尔坐标系
    x = r .* cosd(alpha) .* cosd(theta);
    y = r .* cosd(alpha) .* sind(theta);
    z = r .* sind(alpha);
    
    % 创建拟合对象
    sf = fit([x, y], z, 'poly23');  % 直接使用列向量，不需要转置
    
    % 创建网格点用于绘制曲面
    [xgrid, ygrid] = meshgrid(linspace(min(x), max(x), 100), ...
                             linspace(min(y), max(y), 100));
    zgrid = sf(xgrid, ygrid);
    
    % 创建新图形窗口
    h1 = figure('Position', [100, 100, 800, 600], 'Name', '多项式拟合结果');
    
    % 绘制拟合曲面（添加透明度参数）
    surf(xgrid, ygrid, zgrid, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);  % FaceAlpha=0.5（半透明），EdgeAlpha=0（隐藏边线）
    hold on;
    
    % 添加原始数据点
    scatter3(x, y, z, 50, 'filled', 'r', 'MarkerEdgeColor', 'k');
    
    % 设置图形属性
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title('三维数据点及其多项式拟合曲面');
    grid on;
    colormap('jet');
    colorbar;
    
    % 添加图例
    legend('拟合曲面', '原始数据点', 'Location', 'best');
    
    % 设置视角
    view(-45, 30);
    
    % 获取拟合信息
    [~, gof] = fit([x, y], z, 'poly23');
    
    % 输出拟合统计信息
    fprintf('拟合优度 R-square: %.4f\n', gof.rsquare);
    fprintf('调整后的拟合优度 Adjusted R-square: %.4f\n', gof.adjrsquare);
    fprintf('均方根误差 RMSE: %.4f\n', gof.rmse);
    
    % 输出系数信息
    coeff_names = coeffnames(sf);
    coeff_values = coeffvalues(sf);
    fprintf('\n多项式系数：\n');
    for i = 1:length(coeff_names)
        fprintf('%s: %.4f\n', coeff_names{i}, coeff_values(i));
    end
    
    % 保持图形窗口打开但继续执行
    drawnow;
    pause(1);
    
    %% 第二步：大范围搜索k
    fprintf('\n=== 步骤2：大范围搜索k ===\n');
    
    % 设置k的搜索范围
    k_min = 0.0001;
    k_max = 0.01;
    
    % 使用多个初始值进行优化，避免局部最优
    k_starts = linspace(k_min, k_max, 10);
    best_rmse = inf;
    k_opt_initial = k_min;
    
    % 对每个初始值进行优化
    for k_start = k_starts
        [k_opt, rmse] = fminbnd(@(k) calculate_rmse(k, r, alpha, theta), ...
                               max(k_min, k_start/2), min(k_max, k_start*2));
        if rmse < best_rmse
            best_rmse = rmse;
            k_opt_initial = k_opt;
        end
    end

    % 用初步最优k（k_opt_initial）计算theta值，并存为全局变量（供步骤5调用）
    [theta_wide, valid_points_wide] = calculate_theta(k_opt_initial, r, alpha, theta);
    % 可选：保存对应的RMSE，方便后续对比
    rmse_wide = best_rmse;
    
    fprintf('大范围搜索结果：\n');
    fprintf('初步最优k值: %.6f\n', k_opt_initial);
    fprintf('初步最小RMSE: %.4f°\n', best_rmse);
    
    % 绘制搜索结果
    h2 = figure('Position', [100, 100, 800, 600], 'Name', '大范围k值搜索结果');
    k_range = linspace(k_min, k_max, 100);
    rmse_values = zeros(size(k_range));
    for i = 1:length(k_range)
        rmse_values(i) = calculate_rmse(k_range(i), r, alpha, theta);
    end
    
    plot(k_range, rmse_values, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(k_opt_initial, best_rmse, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    xlabel('k值');
    ylabel('RMSE (度)');
    title('RMSE随k值变化的曲线（大范围搜索）');
    grid on;
    text(k_opt_initial, best_rmse, sprintf('  k = %.6f\n  RMSE = %.4f°', k_opt_initial, best_rmse), ...
        'VerticalAlignment', 'bottom', 'FontSize', 12);
    
    drawnow;
    pause(1);
    
    %% 第三步：解析式拟合曲面
    fprintf('\n=== 步骤3：解析式拟合曲面 ===\n');
    fprintf('使用最优k值: %.6f\n', k_opt_initial);
    
    % 创建更密集的网格点
    alpha_range = linspace(-90, 90, 300);  % 度
    r_range = linspace(0, 400, 300);      % 米
    [alpha_grid, r_grid] = meshgrid(alpha_range, r_range);
    
    % 转换为弧度
    alpha_rad = deg2rad(alpha_grid);
    
    % 初始化两个theta网格，分别存储正根和负根
    theta_grid_pos = nan(size(alpha_grid));
    theta_grid_neg = nan(size(alpha_grid));
    
    % 对每个网格点计算theta
    for i = 1:size(alpha_grid, 1)
        for j = 1:size(alpha_grid, 2)
            kr = k_opt_initial * r_grid(i,j);
            cos_alpha = cos(alpha_rad(i,j));
            sin_alpha = sin(alpha_rad(i,j));
            
            % 计算判别式
            discriminant = 1 - kr * (2*sin_alpha + kr*cos_alpha^2);
            
            if discriminant >= 0 && abs(cos_alpha) > 1e-10
                sqrt_term = sqrt(discriminant);
                denom = kr * cos_alpha;
                
                % 计算两个可能的解
                num1 = 1 + sqrt_term;
                num2 = 1 - sqrt_term;
                
                theta1 = rad2deg(atan2(num1, denom));
                theta2 = rad2deg(atan2(num2, denom));
                
                % 规范化角度
                theta1 = mod(theta1 + 360, 360);
                theta2 = mod(theta2 + 360, 360);
                if theta1 > 90 && theta1 <= 270
                    theta1 = 180 - theta1;
                elseif theta1 > 270
                    theta1 = theta1 - 360;
                end
                if theta2 > 90 && theta2 <= 270
                    theta2 = 180 - theta2;
                elseif theta2 > 270
                    theta2 = theta2 - 360;
                end
                
                % 根据解的大小分配到正负根网格
                if theta1 >= 0 && theta1 <= 90 && theta2 >= 0 && theta2 <= 90
                    if theta1 >= theta2
                        theta_grid_pos(i,j) = theta1;
                        theta_grid_neg(i,j) = theta2;
                    else
                        theta_grid_pos(i,j) = theta2;
                        theta_grid_neg(i,j) = theta1;
                    end
                elseif theta1 >= 0 && theta1 <= 90
                    theta_grid_pos(i,j) = theta1;
                elseif theta2 >= 0 && theta2 <= 90
                    theta_grid_pos(i,j) = theta2;
                end
            end
        end
    end
    
    % 创建figure窗口
    h3 = figure('Position', [100, 100, 800, 600], 'Name', '解析式拟合曲面');
    
    % 绘制两个解的曲面
    h_pos = surf(alpha_grid, r_grid, theta_grid_pos, 'FaceAlpha', 0.7);
    hold on;
    h_neg = surf(alpha_grid, r_grid, theta_grid_neg, 'FaceAlpha', 0.7);
    
    % 设置曲面属性
    set(h_pos, 'EdgeColor', 'none');
    set(h_neg, 'EdgeColor', 'none');
    
    % 使用不同的颜色方案区分两个解
    colormap('jet');
    set(h_pos, 'FaceColor', 'r');  % 正根用红色
    set(h_neg, 'FaceColor', 'b');  % 负根用蓝色
    
    % 添加实验数据点
    scatter3(alpha, r, theta, 100, 'filled', 'g', 'MarkerEdgeColor', 'k');
    
    % 设置图的属性
    xlabel('\alpha (仰角/度)', 'FontSize', 12);
    ylabel('r (距离/m)', 'FontSize', 12);
    zlabel('\theta (方位角/度)', 'FontSize', 12);
    title(['θ(r, α)解析解曲面 (k = ' num2str(k_opt_initial) ')'], 'FontSize', 14);
    view(45, 30);
    grid on;
    
    % 设置坐标轴范围
    xlim([min(alpha_range), max(alpha_range)]);
    ylim([min(r_range), max(r_range)]);
    zlim([0, 90]);
    
    % 添加图例
    legend([h_pos, h_neg], '大角度解', '小角度解', 'Location', 'northeast');
    
    % 优化显示效果
    set(gcf, 'Color', 'white');
    
    % 添加说明文字
    text(min(alpha_range), max(r_range), 90, ...
         {'红色曲面：大角度解', '蓝色曲面：小角度解', '绿色点：实验数据'}, ...
         'FontSize', 10, 'VerticalAlignment', 'top');
    
    drawnow;
    pause(1);
    
    %% 第四步：精细k值搜索
    fprintf('\n=== 步骤4：精细k值搜索 ===\n');
    
    % 在k_opt_initial附近进行精细搜索
    k_center = k_opt_initial;
    k_range = 0.0002;  % 搜索范围缩小到±0.0001
    num_points = 300;  % 增加搜索点数以提高精度
    
    % 创建精细的k值网格
    k_fine = linspace(k_center - k_range/2, k_center + k_range/2, num_points);
    rmse_fine = zeros(size(k_fine));
    theta_calc_best = zeros(size(theta));
    
    % 计算每个k值对应的RMSE
    fprintf('开始精细搜索k值，范围：[%.6f, %.6f]\n', k_fine(1), k_fine(end));
    for i = 1:length(k_fine)
        [rmse_fine(i), theta_calc] = basic_theta_fit(k_fine(i));
        if i == 1 || rmse_fine(i) < min(rmse_fine(1:i-1))
            theta_calc_best = theta_calc;
        end
    end
    
    % 找出最佳k值
    [best_rmse, best_idx] = min(rmse_fine);
    k_opt_final = k_fine(best_idx);
    
    % 输出结果
    fprintf('\n精细搜索结果：\n');
    fprintf('最终最优k值: %.6f\n', k_opt_final);
    fprintf('最终最小RMSE: %.4f°\n', best_rmse);
    fprintf('相对于初始k值的改进：%.2f%%\n', (best_rmse - best_rmse) / best_rmse * 100);
    
    % 绘制精细搜索结果
    h4 = figure('Position', [100, 100, 800, 600], 'Name', '精细k值搜索结果');
    
    % 绘制RMSE曲线
    plot(k_fine, rmse_fine, 'b-', 'LineWidth', 1.5);
    hold on;
    
    % 标记最优点
    plot(k_opt_final, best_rmse, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    % 设置图形属性
    xlabel('k值');
    ylabel('RMSE (度)');
    title('RMSE随k值变化的曲线（精细搜索）');
    grid on;
    
    % 添加文本标注
    text(k_opt_final, best_rmse, sprintf('  k = %.6f\n  RMSE = %.4f°', k_opt_final, best_rmse), ...
        'VerticalAlignment', 'bottom', 'FontSize', 12);
    
    % 设置坐标轴范围，使最优点更容易观察
    k_margin = k_range/10;
    xlim([k_center - k_range/2 - k_margin, k_center + k_range/2 + k_margin]);
    ylim_current = ylim;
    ylim([ylim_current(1) - 0.1, ylim_current(2) + 0.1]);
    
    drawnow;
    pause(1);
    
    %% 第五步：理论曲面和残差分析
    fprintf('\n=== 步骤5：理论曲面和残差分析 ===\n');
    
    % 使用最优k值计算理论theta值
    [theta_theoretical, valid_points] = calculate_theta(k_opt_final, r, alpha, theta);
    
    % 计算残差
    residuals = theta - theta_theoretical;
    valid_residuals = residuals(valid_points);
    
    %保存精细k值搜索后的计算结果
    theta_spic = theta_calc_best;

    % 输出残差统计信息
    fprintf('残差统计：\n');
    fprintf('均值: %.4f°\n', mean(valid_residuals));
    fprintf('标准差: %.4f°\n', std(valid_residuals));
    fprintf('最大残差: %.4f°\n', max(abs(valid_residuals)));
    fprintf('最小残差: %.4f°\n', min(abs(valid_residuals)));
    
    % 创建一个大窗口包含所有图
    h5 = figure('Position', [100, 100, 800, 600], 'Name', '理论曲面拟合与残差分析');
    
    % 创建子图布局
    % 左上：测量值vs计算值
    subplot(2, 2, 1);
    plot(theta, theta_wide, 'b.', 'MarkerSize', 20);
    hold on;
    plot([0, 90], [0, 90], 'r--');
    hold on;
    plot(theta, theta_calc_best, 'r.', 'MarkerSize', 20);
    xlabel('测量角度 (度)');
    ylabel('计算角度 (度)');
    title('(A) 测量角度 vs 计算角度');
    legend('大范围k值搜索计算结果', '理想线', '精细k值搜索计算结果', 'Location', 'best');
    grid on;
    
    % 右上：残差直方图
    subplot(2, 2, 2);
    histogram(valid_residuals, 20, 'Normalization', 'probability');
    xlabel('残差 (度)');
    ylabel('频率');
    title('(B) 残差分布直方图');
    grid on;
    
    % 下方：理论曲面与实验数据对比（占据下方所有空间）
    subplot(2, 2, [3 4]);
    
    % 创建更密集的网格点用于绘制曲面
    alpha_range = linspace(-90, 90, 300);  % 扩大仰角范围，增加密度
    r_range = linspace(0, max(r)*1.2, 300);  % 扩大距离范围，增加密度
    [r_grid, alpha_grid] = meshgrid(r_range, alpha_range);
    
    % 初始化theta网格
    theta_grid_pos = nan(size(r_grid));
    theta_grid_neg = nan(size(r_grid));
    
    % 对每个网格点计算theta
    for i = 1:size(r_grid, 1)
        for j = 1:size(r_grid, 2)
            kr = k_opt_final * r_grid(i,j);
            alpha_rad = deg2rad(alpha_grid(i,j));
            cos_alpha = cos(alpha_rad);
            sin_alpha = sin(alpha_rad);
            
            % 计算判别式
            discriminant = 1 - kr * (2*sin_alpha + kr*cos_alpha^2);
            
            if discriminant >= 0 && abs(cos_alpha) > 1e-10
                sqrt_term = sqrt(discriminant);
                denom = kr * cos_alpha;
                
                % 计算两个可能的解
                num1 = 1 + sqrt_term;
                num2 = 1 - sqrt_term;
                
                theta1 = rad2deg(atan2(num1, denom));
                theta2 = rad2deg(atan2(num2, denom));
                
                % 规范化角度
                theta1 = mod(theta1 + 360, 360);
                theta2 = mod(theta2 + 360, 360);
                if theta1 > 90 && theta1 <= 270
                    theta1 = 180 - theta1;
                elseif theta1 > 270
                    theta1 = theta1 - 360;
                end
                if theta2 > 90 && theta2 <= 270
                    theta2 = 180 - theta2;
                elseif theta2 > 270
                    theta2 = theta2 - 360;
                end
                
                % 根据解的大小分配到正负根网格
                if theta1 >= 0 && theta1 <= 90 && theta2 >= 0 && theta2 <= 90
                    if theta1 >= theta2
                        theta_grid_pos(i,j) = theta1;
                        theta_grid_neg(i,j) = theta2;
                    else
                        theta_grid_pos(i,j) = theta2;
                        theta_grid_neg(i,j) = theta1;
                    end
                elseif theta1 >= 0 && theta1 <= 90
                    theta_grid_pos(i,j) = theta1;
                elseif theta2 >= 0 && theta2 <= 90
                    theta_grid_pos(i,j) = theta2;
                end
            end
        end
    end
    
    % 绘制两个解的曲面
    h_pos = surf(r_grid, alpha_grid, theta_grid_pos, 'FaceAlpha', 0.7);
    hold on;
    h_neg = surf(r_grid, alpha_grid, theta_grid_neg, 'FaceAlpha', 0.7);
    
    % 设置曲面属性
    set(h_pos, 'EdgeColor', 'none');
    set(h_neg, 'EdgeColor', 'none');
    
    % 使用不同的颜色方案区分两个解
    colormap('jet');
    set(h_pos, 'FaceColor', 'r');  % 正根用红色
    set(h_neg, 'FaceColor', 'b');  % 负根用蓝色
    
    % 添加实验数据点
    scatter3(r, alpha, theta, 100, 'filled', 'g', 'MarkerEdgeColor', 'k');
    
    % 设置图形属性
    xlabel('距离 r (m)');
    ylabel('仰角 \alpha (度)');
    zlabel('方位角 \theta (度)');
    title(['(C) 理论曲面与实验数据对比 (k = ' num2str(k_opt_final) ')']);
    grid on;
    
    % 设置坐标轴范围
    xlim([0, max(r)*1.2]);
    ylim([-90, 90]);
    zlim([0, 90]);
    
    % 添加图例
    legend([h_pos, h_neg], '大角度解', '小角度解', 'Location', 'northeast');
    
    % 优化视角
    view(45, 30);
    
    % 优化显示效果
    set(gcf, 'Color', 'white');
    
    % 调整子图间距
    set(gcf, 'Units', 'normalized');
    set(findall(gcf, 'Type', 'axes'), 'FontSize', 10);
    
    drawnow;
    pause(1);
    
    %% 第六步：k函数优化
    fprintf('\n=== 步骤6：k函数优化 ===\n');
    
    % 使用之前找到的最优k值作为基准
    k_base = k_opt_final;
    fprintf('基准k值: %.6f\n', k_base);
    
    % 初始化扰动参数 [k_base, a, b, c, d, e]
    % k = k_base * (1 + a*r + b*alpha + c*r*alpha + d*r^2 + e*alpha^2)
    params = [k_base, 0, 0, 0, 0, 0];  % [k_base, a, b, c, d, e]
    
    % 设置优化选项
    options = optimset('Display', 'iter', ...
                      'MaxFunEvals', 3000, ...
                      'MaxIter', 1500, ...
                      'TolX', 1e-8, ...
                      'TolFun', 1e-8);
    
    % 定义参数边界
    lb = [0.002, -0.001, -0.001, -0.0001, -0.0001, -0.0001];  % 下界
    ub = [0.003,  0.001,  0.001,  0.0001,  0.0001,  0.0001];  % 上界
    
    % 使用模式搜索算法进行优化
    fprintf('开始优化扰动参数...\n');
    [optimal_params, fval] = patternsearch(@(p) calculate_rmse_with_perturbation(p, r, alpha, theta), ...
                                         params, [], [], [], [], lb, ub, [], options);
    
    % 使用fminsearch进行局部优化
    [optimal_params, fval] = fminsearch(@(p) calculate_rmse_with_perturbation(p, r, alpha, theta), ...
                                      optimal_params, options);
    
    % 计算最终结果
    [final_rmse, theta_calc, k_values] = calculate_rmse_with_perturbation(optimal_params, r, alpha, theta);
    
    % 输出优化结果
    fprintf('\n优化结果：\n');
    fprintf('k_base = %.6f\n', optimal_params(1));
    fprintf('a (r系数) = %.6f\n', optimal_params(2));
    fprintf('b (alpha系数) = %.6f\n', optimal_params(3));
    fprintf('c (r*alpha系数) = %.6f\n', optimal_params(4));
    fprintf('d (r^2系数) = %.6f\n', optimal_params(5));
    fprintf('e (alpha^2系数) = %.6f\n', optimal_params(6));
    fprintf('RMSE = %.4f°\n', final_rmse);
    
    % 计算误差统计
    errors = theta_calc - theta;
    valid_errors = errors(~isnan(errors));
    fprintf('\n误差统计：\n');
    fprintf('均值: %.4f°\n', mean(valid_errors));
    fprintf('标准差: %.4f°\n', std(valid_errors));
    fprintf('最大误差: %.4f°\n', max(abs(valid_errors)));
    fprintf('最小误差: %.4f°\n', min(abs(valid_errors)));

    %保存k函数优化的计算结果
    theta_kfun = theta_calc;
    
    % 创建figure显示优化结果
    h6 = figure('Position', [100, 100, 800, 600], 'Name', 'k函数优化结果');
    
    % 子图1：测量值vs计算值
    subplot(2, 2, 1);
    plot(theta, theta_calc, 'r.', 'MarkerSize', 20);
    hold on;
    plot([0, 90], [0, 90], 'r--');
    hold on;
    plot(theta, theta_spic, 'b.', 'MarkerSize', 20);
    xlabel('测量角度 (度)');
    ylabel('计算角度 (度)');
    title('(A) 测量角度 vs 计算角度');
    legend('基于k函数的计算结果', '理想线', '精细k值搜索计算结果', 'Location', 'best');
    grid on;
    
    % 子图2：残差分布
    subplot(2, 2, 2);
    histogram(valid_errors, 20, 'Normalization', 'probability');
    xlabel('残差 (度)');
    ylabel('频率');
    title(sprintf('(B) 残差分布 (RMSE=%.4f°)', final_rmse));
    grid on;
    
    % 子图3：k值与距离的关系
    subplot(2, 2, 3);
    scatter(r, k_values, 50, theta, 'filled');
    colorbar;
    xlabel('距离 r (m)');
    ylabel('k值');
    title('(C) k值与距离的关系 (颜色表示theta)');
    grid on;
    
    % 子图4：k值与alpha的关系
    subplot(2, 2, 4);
    scatter(alpha, k_values, 50, theta, 'filled');
    colorbar;
    xlabel('仰角 alpha (度)');
    ylabel('k值');
    title('(D) k值与仰角的关系 (颜色表示theta)');
    grid on;
    
    % 调整布局
    set(gcf, 'Color', 'white');
    set(findall(gcf, 'Type', 'axes'), 'FontSize', 10);
    
    % 保存第六步的优化结果
    step6_results = struct();
    step6_results.optimal_params = optimal_params;
    step6_results.final_rmse = final_rmse;
    step6_results.theta_calc = theta_calc;
    step6_results.k_values = k_values;
    save('step6_results.mat', '-struct', 'step6_results');
    
    drawnow;
    pause(1);
    
    %% 第七步：theta表达式扰动优化
    fprintf('\n=== 步骤7：theta表达式扰动优化 ===\n');
    
    % 加载第六步的优化结果
    step6_results = load('step6_results.mat');
    k_params = step6_results.optimal_params;
    
    % 输出k函数参数
    fprintf('使用的k函数参数：\n');
    fprintf('k_base = %.6f\n', k_params(1));
    fprintf('a (r系数) = %.6f\n', k_params(2));
    fprintf('b (alpha系数) = %.6f\n', k_params(3));
    fprintf('c (r*alpha系数) = %.6f\n', k_params(4));
    fprintf('d (r^2系数) = %.6f\n', k_params(5));
    fprintf('e (alpha^2系数) = %.6f\n\n', k_params(6));
    
    % 初始化theta扰动参数 [theta_a1, theta_a2, theta_a3, theta_a4, theta_a5]
    % theta_new = theta + theta_a1 + theta_a2*r + theta_a3*sin(alpha) + 
    %            theta_a4*sin(2*alpha) + theta_a5*r.*sin(alpha)
    theta_params = [0, 0, 0, 0, 0];
    
    % 设置优化选项
    options = optimset('Display', 'iter', ...
                      'MaxFunEvals', 5000, ...
                      'MaxIter', 2500, ...
                      'TolX', 1e-6, ...
                      'TolFun', 1e-6);
    
    % 使用fminsearch进行优化
    fprintf('开始优化theta扰动参数...\n');
    [optimal_theta_params, fval] = fminsearch(@(p) calculate_rmse_with_theta_perturbation(p, k_params, r, alpha, theta), ...
                                            theta_params, options);
    
    % 计算最终结果
    [final_rmse, theta_calc, theta_perturbation] = calculate_rmse_with_theta_perturbation(optimal_theta_params, k_params, r, alpha, theta);
    
    % 输出优化结果
    fprintf('\n优化结果：\n');
    fprintf('theta_a1 (角度扰动) = %.6f\n', optimal_theta_params(1));
    fprintf('theta_a2 (距离扰动) = %.6f\n', optimal_theta_params(2));
    fprintf('theta_a3 (alpha角度扰动) = %.6f\n', optimal_theta_params(3));
    fprintf('theta_a4 (alpha二倍角扰动) = %.6f\n', optimal_theta_params(4));
    fprintf('theta_a5 (r-alpha交互扰动) = %.6f\n', optimal_theta_params(5));
    fprintf('最终RMSE = %.4f°\n', final_rmse);
    
    % 创建figure显示优化结果
    h7 = figure('Position', [100, 100, 800, 600], 'Name', 'theta扰动优化结果');
    
    % 子图1：测量值vs计算值
    subplot(2, 2, 1);
    plot(theta, theta_calc, 'r.', 'MarkerSize', 20);
    hold on;
    plot([0, 90], [0, 90], 'r--');
    hold on;
    plot(theta, theta_kfun, 'b.', 'MarkerSize', 20);
    xlabel('测量角度 (度)');
    ylabel('计算角度 (度)');
    legend('添加theta扰动的计算结果', '理想线', '基于k函数的计算结果', 'Location', 'best');
    xlabel('测量角度 (度)');
    ylabel('计算角度 (度)');
    title('(A) 测量角度 vs 计算角度');
    grid on;
    
    % 子图2：残差分布
    subplot(2, 2, 2);
    errors = theta_calc - theta;
    histogram(errors, 20, 'Normalization', 'probability');
    xlabel('残差 (度)');
    ylabel('频率');
    title(sprintf('(B) 残差分布 (RMSE=%.4f°)', final_rmse));
    grid on;
    
    % 子图3：扰动量与距离的关系
    subplot(2, 2, 3);
    scatter(r, theta_perturbation, 40, alpha, 'filled');
    colorbar;
    xlabel('距离 r (m)');
    ylabel('theta扰动量 (度)');
    title('(C) 扰动量与距离的关系 (颜色表示alpha)');
    grid on;
    
    % 子图4：扰动量与alpha的关系
    subplot(2, 2, 4);
    scatter(alpha, theta_perturbation, 40, r, 'filled');
    colorbar;
    xlabel('仰角 alpha (度)');
    ylabel('theta扰动量 (度)');
    title('(D) 扰动量与仰角的关系 (颜色表示r)');
    grid on;
    
    % 保存优化结果
    step7_results = struct();
    step7_results.optimal_theta_params = optimal_theta_params;
    step7_results.final_rmse = final_rmse;
    step7_results.theta_calc = theta_calc;
    step7_results.theta_perturbation = theta_perturbation;
    save('step7_results.mat', '-struct', 'step7_results');
    
    drawnow;
    pause(1);
    
    %% 第八步：theta校正函数优化
    fprintf('\n=== 步骤8：theta校正函数优化 ===\n');
    
    % 设定分界点
    theta_threshold = 45;
    
    % 分离高低角度数据
    low_angle_idx = theta < theta_threshold;
    high_angle_idx = theta >= theta_threshold;
    
    % 计算需要校正的差值
    delta_theta = theta - theta_calc;
    
    % 创建拟合数据矩阵
    % 低角度使用增强的模型
    X_low = @(theta, r, alpha) [ones(size(theta)), ...
                               theta, theta.^2, ...
                               r/300, (r/300).^2, ...
                               sin(deg2rad(alpha)), cos(deg2rad(alpha)), ...
                               sin(deg2rad(2*alpha)), ...
                               (r/300) .* sin(deg2rad(alpha)), ...
                               alpha/90];  % 添加归一化的alpha线性项
    
    % 高角度使用简化模型
    X_high = @(theta, r, alpha) [ones(size(theta)), theta, theta.^2, ...
                                r/300, (r/300).^2, ...
                                sin(deg2rad(alpha)), sin(deg2rad(2*alpha))];
    
    % 构建设计矩阵
    X_low_mat = X_low(theta_calc(low_angle_idx), r(low_angle_idx), alpha(low_angle_idx));
    X_high_mat = X_high(theta_calc(high_angle_idx), r(high_angle_idx), alpha(high_angle_idx));
    
    % 添加Tikhonov正则化（岭回归）
    lambda_low = 0.005;
    lambda_high = 0.005;
    
    % 使用正则化最小二乘法拟合参数
    I_low = eye(size(X_low_mat, 2));
    I_high = eye(size(X_high_mat, 2));
    
    correction_params_low = (X_low_mat' * X_low_mat + lambda_low * I_low) \ ...
                          (X_low_mat' * delta_theta(low_angle_idx));
    
    correction_params_high = (X_high_mat' * X_high_mat + lambda_high * I_high) \ ...
                           (X_high_mat' * delta_theta(high_angle_idx));
    
    % 计算校正后的theta值
    theta_corrected = theta_calc;
    theta_corrected(low_angle_idx) = theta_calc(low_angle_idx) + ...
                                    X_low_mat * correction_params_low;
    theta_corrected(high_angle_idx) = theta_calc(high_angle_idx) + ...
                                     X_high_mat * correction_params_high;
    
    % 计算RMSE
    rmse = sqrt(mean((theta - theta_corrected).^2));
    rmse_low = sqrt(mean((theta(low_angle_idx) - theta_corrected(low_angle_idx)).^2));
    rmse_high = sqrt(mean((theta(high_angle_idx) - theta_corrected(high_angle_idx)).^2));
    
    % 创建figure显示优化结果
    h8 = figure('Position', [100, 100, 800, 600], 'Name', 'theta校正函数优化结果');
    
    % % 子图1：原始误差和拟合校正
    % subplot(2, 2, 1);
    % plot(theta_calc(low_angle_idx), delta_theta(low_angle_idx), 'b.', 'MarkerSize', 12);
    % hold on;
    % plot(theta_calc(high_angle_idx), delta_theta(high_angle_idx), 'r.', 'MarkerSize', 12);
    % plot(theta_calc(low_angle_idx), X_low_mat * correction_params_low, 'c.', 'MarkerSize', 8);
    % plot(theta_calc(high_angle_idx), X_high_mat * correction_params_high, 'm.', 'MarkerSize', 8);
    % xlabel('计算得到的\theta (度)');
    % ylabel('\Delta\theta (度)');
    % title('(A) 原始误差和拟合校正');
    % legend('低角度误差', '高角度误差', '低角度校正', '高角度校正', 'Location', 'best');
    % grid on;
    % 
    % % 子图2：距离相关性
    % subplot(2, 2, 2);
    % plot(r(low_angle_idx), delta_theta(low_angle_idx), 'b.', 'MarkerSize', 12);
    % hold on;
    % plot(r(high_angle_idx), delta_theta(high_angle_idx), 'r.', 'MarkerSize', 12);
    % xlabel('距离 r (m)');
    % ylabel('\Delta\theta (度)');
    % title('(B) 误差与距离的关系');
    % legend('低角度', '高角度', 'Location', 'best');
    % grid on;
    % 
    % % 子图3：alpha角度相关性
    % subplot(2, 2, 3);
    % plot(alpha(low_angle_idx), delta_theta(low_angle_idx), 'b.', 'MarkerSize', 12);
    % hold on;
    % plot(alpha(high_angle_idx), delta_theta(high_angle_idx), 'r.', 'MarkerSize', 12);
    % xlabel('\alpha (度)');
    % ylabel('\Delta\theta (度)');
    % title('(C) 误差与\alpha的关系');
    % legend('低角度', '高角度', 'Location', 'best');
    % grid on;
    
    % 子图4：校正前后对比
    % subplot(2, 2, 4);
    plot(theta, theta_calc, 'b.', 'MarkerSize', 20);
    hold on;
    plot(theta, theta_corrected, 'r.', 'MarkerSize', 20);
    plot([0 90], [0 90], 'k--');
    xlabel('测量\theta (度)');
    ylabel('计算\theta (度)');
    title('校正前后对比');
    legend('校正前', '校正后', '理想线', 'Location', 'best');
    grid on;
    
    % 输出结果
    fprintf('\n低角度(≤%.1f°)校正函数参数：\n', theta_threshold);
    fprintf('常数项: %.6f\n', correction_params_low(1));
    fprintf('theta项: %.6f\n', correction_params_low(2));
    fprintf('theta^2项: %.6f\n', correction_params_low(3));
    fprintf('r/300项: %.6f\n', correction_params_low(4));
    fprintf('(r/300)^2项: %.6f\n', correction_params_low(5));
    fprintf('sin(alpha)项: %.6f\n', correction_params_low(6));
    fprintf('cos(alpha)项: %.6f\n', correction_params_low(7));
    fprintf('sin(2*alpha)项: %.6f\n', correction_params_low(8));
    fprintf('r*sin(alpha)交互项: %.6f\n', correction_params_low(9));
    fprintf('alpha/90项: %.6f\n', correction_params_low(10));
    
    fprintf('\n高角度(>%.1f°)校正函数参数：\n', theta_threshold);
    fprintf('常数项: %.6f\n', correction_params_high(1));
    fprintf('theta项: %.6f\n', correction_params_high(2));
    fprintf('theta^2项: %.6f\n', correction_params_high(3));
    fprintf('r/300项: %.6f\n', correction_params_high(4));
    fprintf('(r/300)^2项: %.6f\n', correction_params_high(5));
    fprintf('sin(alpha)项: %.6f\n', correction_params_high(6));
    fprintf('sin(2*alpha)项: %.6f\n', correction_params_high(7));
    
    fprintf('\nRMSE改进：\n');
    fprintf('校正前RMSE: %.4f°\n', sqrt(mean((theta - theta_calc).^2)));
    fprintf('校正后总体RMSE: %.4f°\n', rmse);
    fprintf('低角度RMSE: %.4f°\n', rmse_low);
    fprintf('高角度RMSE: %.4f°\n', rmse_high);
    
    % 保存优化结果
    step8_results = struct();
    step8_results.correction_params_low = correction_params_low;
    step8_results.correction_params_high = correction_params_high;
    step8_results.theta_threshold = theta_threshold;
    step8_results.rmse = rmse;
    step8_results.rmse_low = rmse_low;
    step8_results.rmse_high = rmse_high;
    step8_results.theta_corrected = theta_corrected;
    save('step8_results.mat', '-struct', 'step8_results');
    
    % 创建校正函数文件
    create_correction_function(correction_params_low, correction_params_high, theta_threshold);
    
    drawnow;
    pause(1);
end

function create_correction_function(params_low, params_high, theta_threshold)
    % 创建校正函数文件
    filename = 'theta_correction_function.m';
    fid = fopen(filename, 'w');
    
    fprintf(fid, 'function delta = theta_correction_function(theta, r, alpha)\n');
    fprintf(fid, '    %% theta校正函数\n');
    fprintf(fid, '    %% theta: 计算得到的theta值（度）\n');
    fprintf(fid, '    %% r: 距离（米）\n');
    fprintf(fid, '    %% alpha: 角度（度）\n');
    fprintf(fid, '    \n');
    fprintf(fid, '    %% 判断使用哪组参数\n');
    fprintf(fid, '    use_low = theta <= %.1f;\n', theta_threshold);
    fprintf(fid, '    \n');
    fprintf(fid, '    %% 初始化输出\n');
    fprintf(fid, '    delta = zeros(size(theta));\n');
    fprintf(fid, '    \n');
    fprintf(fid, '    %% 低角度校正（增强模型）\n');
    fprintf(fid, '    if any(use_low)\n');
    fprintf(fid, '        delta(use_low) = %.6f + ...\n', params_low(1));
    fprintf(fid, '            %.6f * theta(use_low) + ...\n', params_low(2));
    fprintf(fid, '            %.6f * theta(use_low).^2 + ...\n', params_low(3));
    fprintf(fid, '            %.6f * r(use_low)/300 + ...\n', params_low(4));
    fprintf(fid, '            %.6f * (r(use_low)/300).^2 + ...\n', params_low(5));
    fprintf(fid, '            %.6f * sin(deg2rad(alpha(use_low))) + ...\n', params_low(6));
    fprintf(fid, '            %.6f * cos(deg2rad(alpha(use_low))) + ...\n', params_low(7));
    fprintf(fid, '            %.6f * sin(deg2rad(2*alpha(use_low))) + ...\n', params_low(8));
    fprintf(fid, '            %.6f * (r(use_low)/300) .* sin(deg2rad(alpha(use_low))) + ...\n', params_low(9));
    fprintf(fid, '            %.6f * alpha(use_low)/90;\n', params_low(10));
    fprintf(fid, '    end\n');
    fprintf(fid, '    \n');
    fprintf(fid, '    %% 高角度校正（简化模型）\n');
    fprintf(fid, '    if any(~use_low)\n');
    fprintf(fid, '        delta(~use_low) = %.6f + ...\n', params_high(1));
    fprintf(fid, '            %.6f * theta(~use_low) + ...\n', params_high(2));
    fprintf(fid, '            %.6f * theta(~use_low).^2 + ...\n', params_high(3));
    fprintf(fid, '            %.6f * r(~use_low)/300 + ...\n', params_high(4));
    fprintf(fid, '            %.6f * (r(~use_low)/300).^2 + ...\n', params_high(5));
    fprintf(fid, '            %.6f * sin(deg2rad(alpha(~use_low))) + ...\n', params_high(6));
    fprintf(fid, '            %.6f * sin(deg2rad(2*alpha(~use_low)));\n', params_high(7));
    fprintf(fid, '    end\n');
    fprintf(fid, 'end\n');
    
    fclose(fid);
end

function rmse = calculate_rmse(k, r, alpha, theta_measured)  %通过给定的k用解析式计算出theta并减去测量值计算rmse
    % 使用向量化计算theta
    [theta_predicted, ~] = calculate_theta(k, r, alpha, theta_measured);
    
    % 计算RMSE
    residuals = theta_measured - theta_predicted;
    valid_residuals = residuals(~isnan(residuals));
    rmse = sqrt(mean(valid_residuals.^2));
end

function [rmse, theta_calc] = basic_theta_fit(k_value)   %详细版的针对给定k计算rmse的函数，比上面的多输出了每个k对应的计算结果
    % 使用全局数据
    global DATA;
    r = DATA(:,1);
    alpha = DATA(:,2);
    theta_measured = DATA(:,3);
    
    % 使用向量化计算theta
    [theta_calc, valid_points] = calculate_theta(k_value, r, alpha, theta_measured);
    
    % 计算RMSE
    rmse = sqrt(mean((theta_measured(valid_points) - theta_calc(valid_points)).^2));
    
    % 输出详细信息
    fprintf('\n计算结果：\n');
    fprintf('k值: %.6f\n', k_value);
    fprintf('RMSE: %.4f°\n', rmse);
    fprintf('\n详细对比：\n');
    fprintf('序号\t   r(m)\t alpha(°)\t theta测量(°)\t theta计算(°)\t 误差(°)\n');
    for i = 1:length(r)
        fprintf('%3d\t%7.2f\t%8.2f\t%12.2f\t%12.2f\t%8.2f\n', ...
                i, r(i), alpha(i), theta_measured(i), theta_calc(i), ...
                theta_calc(i) - theta_measured(i));
    end
end

function [theta_calc, valid_points] = calculate_theta(k, r, alpha, theta_measured)    %根据给定的k,r,alpha计算theta预测值，再结合实际测量值给出更接近实测的theta值
    % 向量化计算theta
    alpha_rad = deg2rad(alpha);
    cos_alpha = cos(alpha_rad);
    sin_alpha = sin(alpha_rad);
    kr = k * r;
    
    % 计算判别式
    discriminant = 1 - kr .* (2*sin_alpha + kr.*cos_alpha.^2);
    
    % 初始化输出
    theta_calc = zeros(size(theta_measured));
    valid_points = discriminant >= 0 & abs(cos_alpha) > 1e-10;
    
    % 只对有效点计算theta
    if any(valid_points)
        sqrt_term = sqrt(discriminant(valid_points));
        denom = kr(valid_points) .* cos_alpha(valid_points);
        
        % 计算两个可能的解
        theta1 = rad2deg(atan((1 + sqrt_term) ./ denom));
        theta2 = rad2deg(atan((1 - sqrt_term) ./ denom));
        
        % 选择更接近测量值的解
        diff1 = abs(theta1 - theta_measured(valid_points));
        diff2 = abs(theta2 - theta_measured(valid_points));
        
        % 初始化为theta1
        theta_calc(valid_points) = theta1;
        
        % 找出theta2更好的点
        use_theta2 = diff2 < diff1;
        valid_indices = find(valid_points);
        theta_calc(valid_indices(use_theta2)) = theta2(use_theta2);
    end
    
    % 将无效点设为NaN
    theta_calc(~valid_points) = nan;
end


function [rmse, theta_calc, k_values] = calculate_rmse_with_perturbation(params, r, alpha, theta)    %基于k扰动模型计算theta和对应的RMSE，输出每个r和alpha对应的最小误差k值
    % 提取参数
    k_base = params(1);
    a = params(2);
    b = params(3);
    c = params(4);
    d = params(5);
    e = params(6);
    
    % 计算每个点的k值，添加交叉项和二次项
    k_values = k_base * (1 + a*r + b*alpha + c*r.*alpha + d*r.^2 + e*alpha.^2);
    
    % 初始化
    theta_calc = zeros(size(theta));
    alpha_rad = deg2rad(alpha);
    
    % 计算每个点的theta
    for i = 1:length(r)
        % 计算判别式
        kr = k_values(i) * r(i);
        cos_alpha = cos(alpha_rad(i));
        sin_alpha = sin(alpha_rad(i));
        discriminant = 1 - kr * (2*sin_alpha + kr*cos_alpha^2);
        
        if discriminant >= 0 && abs(cos_alpha) > 1e-10
            % 计算两个可能的解
            sqrt_term = sqrt(discriminant);
            denom = kr * cos_alpha;
            
            theta1 = rad2deg(atan((1 + sqrt_term) / denom));
            theta2 = rad2deg(atan((1 - sqrt_term) / denom));
            
            % 选择更接近测量值的解
            if abs(theta1 - theta(i)) < abs(theta2 - theta(i))
                theta_calc(i) = theta1;
            else
                theta_calc(i) = theta2;
            end
        else
            % 无物理解时，使用测量值
            theta_calc(i) = theta(i);
        end
    end
    
    % 计算RMSE
    rmse = sqrt(mean((theta - theta_calc).^2));
end

function [rmse, theta_calc, theta_perturbation] = calculate_rmse_with_theta_perturbation(theta_params, k_params, r, alpha, theta)   %基于k扰动模型和theta扰动模型计算RMSE
    % 提取theta扰动参数
    theta_a1 = theta_params(1);  % 角度扰动
    theta_a2 = theta_params(2);  % 距离扰动
    theta_a3 = theta_params(3);  % alpha角度扰动
    theta_a4 = theta_params(4);  % alpha二倍角扰动
    theta_a5 = theta_params(5);  % r-alpha交互扰动
    
    % 提取k函数参数
    k_base = k_params(1);
    a = k_params(2);
    b = k_params(3);
    c = k_params(4);
    d = k_params(5);
    e = k_params(6);
    
    % 计算k值
    k_values = k_base * (1 + a*r + b*alpha + c*r.*alpha + d*r.^2 + e*alpha.^2);
    
    % 初始化
    theta_calc = zeros(size(theta));
    alpha_rad = deg2rad(alpha);
    
    % 计算每个点的theta
    for i = 1:length(r)
        % 计算判别式
        kr = k_values(i) * r(i);
        cos_alpha = cos(alpha_rad(i));
        sin_alpha = sin(alpha_rad(i));
        discriminant = 1 - kr * (2*sin_alpha + kr*cos_alpha^2);
        
        if discriminant >= 0 && abs(cos_alpha) > 1e-10
            % 计算两个可能的解
            sqrt_term = sqrt(discriminant);
            denom = kr * cos_alpha;
            
            theta1 = rad2deg(atan((1 + sqrt_term) / denom));
            theta2 = rad2deg(atan((1 - sqrt_term) / denom));
            
            % 选择更接近测量值的解
            if abs(theta1 - theta(i)) < abs(theta2 - theta(i))
                theta_calc(i) = theta1;
            else
                theta_calc(i) = theta2;
            end
        else
            % 无物理解时，使用测量值
            theta_calc(i) = theta(i);
        end
    end
    
    % 添加扰动项
    theta_perturbation = theta_a1 + theta_a2*r + theta_a3*sin(deg2rad(alpha)) + ...
                        theta_a4*sin(deg2rad(2*alpha)) + theta_a5*r.*sin(deg2rad(alpha));
    theta_calc = theta_calc + theta_perturbation;
    
    % 计算RMSE
    rmse = sqrt(mean((theta - theta_calc).^2));
end
