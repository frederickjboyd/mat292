function [t, y] = adaptiveEuler(f, t_0, t_N, y_0, h)
    % Implementation of the Improved Euler Method
    % f: Function
    % t_0: Start point
    % t_N: End point
    % y_0: Initial condition
    % h: Step size
    
    % Setting parameters
    tol = 1e-8;
    
    % Setting up initial condition
    t_current = t_0;
    y_current = y_0;
    t = t_0;
    y = y_0;
    
    % Storing initial step size
    h_temp = h;
    
    while (t_current < t_N)
        [D, Z] = helper(f, t_current, y_current, h_temp);
        L = length(y);
        t_current = t(L);
        y_current = y(L);
        
        if (abs(D) < tol)
            t_current = t_current + h_temp;
            y_current = Z + D;
            t = horzcat(t, t_current);
            y = horzcat(y, y_current);
            h_temp = h;
        else
            h_temp = 0.9 * h_temp * min(max(tol / abs(D), 0.3), 2);
        end
    end
end

function [D, Z] = helper(f, t, y, h)
    % Finding slope of current point
    k_1 = f(t, y);
    
    % Calculating position of next point in 1 step
    Y = y + (h * k_1);
    
    % Calculating position of next point in 2 half-steps
    Z_1 = y + ((h/2) * k_1);
    k_2 = f(t + (h/2), Z_1);
    Z = Z_1 + ((h/2) * k_2);

    % Calculating difference between 1 step and 2 step methods
    D = Z - Y;
end
