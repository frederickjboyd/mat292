function [t, y] = DE2_boydfred(p, q, g, t0, tN, y0, y1, h)
    % y'' + p(t)*y' + q(t)*y = g(t)
    % y'' = -p(t)*y' - q(t)*y + g(t)

    % Calculating number of points to use
    N = ((tN - t0) / h) + 1;
    N = round(N);
    
    % Allocating space for t, y vectors
    t = linspace(t0, tN, N);
    y = zeros(1, N);
    
    % Setting initial conditions
    y(1) = y0;
    k_prev = y1; % Storing previous slope
    
    % Calculating first step using given slope
    y(2) = y(1) + y1 * h;
    
    % Loop
    for i = 2:N
        % Current t value
        t_curr = t(i);
        % Previous y value
        y_prev = y(i-1);
        % Value of y''
        k_slope = -p(t_curr) * k_prev - q(t_curr) * y_prev + g(t_curr);
        % Value of y'
        k = k_prev + k_slope * h;
        % Storing current slope value in order to calculate slope for next loop
        k_prev = k;
        % Calculating next y value
        y(i) = y_prev + k * h;
    end
end
