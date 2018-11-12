function [t, x] = solvesystem_boydfred(f_1, f_2, t_0, t_N, x_0, h)
    % Calculating number of points to use
    N = ((t_N - t_0) / h) + 1;
    N = round(N);
    
    % Allocating space for t, y vectors
    t = linspace(t_0, t_N, N);
    x = zeros(2, N);
    
    % Setting initial condition
    x(1, 1) = x_0(1);
    x(2, 1) = x_0(2);
    
    for i = 2:N
        % Finding next slope using slope of current points
        % k_1: slopes for x_1
        % k_2: slopes for x_2
        k_1_1 = f_1(x(1, i - 1), x(2, i - 1));
        k_2_1 = f_2(x(1, i - 1), x(2, i - 1));
        u_1_next = x(1, i - 1) + (h * k_1_1);
        u_2_next = x(2, i - 1) + (h * k_2_1);
        % Finding slope of next points
        k_1_2 = f_1(u_1_next, u_2_next);
        k_2_2 = f_2(u_1_next, u_2_next);
        % Averaging two slopes
        k_1_avg = (k_1_1 + k_1_2) / 2;
        k_2_avg = (k_2_1 + k_2_2) / 2;
        % Calculating next point
        x(1, i) = x(1, i - 1) + (h * k_1_avg);
        x(2, i) = x(2, i - 1) + (h * k_2_avg);
    end
end
