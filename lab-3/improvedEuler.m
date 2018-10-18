function [t, y] = improvedEuler(f, t_0, t_N, y_0, h)
    % Implementation of the Improved Euler Method
    % f: Function
    % t_0: Start point
    % t_N: End point
    % y_0: Initial condition
    % h: Step size
    
    % Calculating number of points to use
    N = ((t_N - t_0) / h) + 1;
    N = round(N);
    
    % Allocating space for t, y vectors
    t = linspace(t_0, t_N, N);
    y = zeros(1, N);
    
    % Setting initial condition
    y(1) = y_0;
    
    for i = 2:N
        % Finding next point using slope of previous point
        k_1 = f(t(i - 1), y(i - 1));
        u_next = y(i - 1) + (h * k_1);
        % Finding slope of next point and averaging it with first slope
        k_2 = f(t(i), u_next);
        k_avg = (k_1 + k_2) / 2;
        % Calculating next point
        y(i) = y(i - 1) + (h * k_avg);
    end
end
