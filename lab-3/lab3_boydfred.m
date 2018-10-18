%% ODE Lab: Creating your own ODE solver in MATLAB
%
% In this lab, you will write your own ODE solver for the Improved Euler 
% method (also known as the Heun method), and compare its results to those 
% of |ode45|.
%
% You will also learn how to write a function in a separate m-file and 
% execute it.
% 
% Opening the m-file lab3.m in the MATLAB editor, step through each
% part using cell mode to see the results.  Compare the output with the
% PDF, which was generated from this m-file.
%
% There are six (6) exercises in this lab that are to be handed in on the
% due date. Write your solutions in the template, including
% appropriate descriptions in each step. Save the .m files and submit them 
% online on Quercus.
%
% MAT292, Fall 2018, Stinchcombe & Khovanskii, modified from
% MAT292, Fall 2017, Stinchcombe & Sinnamon, modified from
% MAT292, Fall 2015, Sousa, modified from
% MAT292, Fall 2013, Sinnamon & Sousa, modified from
% MAT292, Fall 2011, Hart & Pym

%% Student Information
%
% Student Name: Frederick Boyd
%
% Student Number: 
%

%% Creating new functions using m-files.
%  
% Create a new function in a separate m-file:
%
% Specifics:  Create a text file with the file name f.m
% with the following lines of code (text):
%
%  function y = f(a,b,c) 
%  y = a+b+c;
%
% Now MATLAB can call the new function f (which simply accepts 3 numbers
% and adds them together).  
% To see how this works, type the following in the matlab command window:
% sum = f(1,2,3)

%% Exercise 1
%
% Objective: Write your own ODE solver (using the Heun/Improved Euler
% Method).
%
% Details: This m-file should be a function which accepts as variables 
% (t0,tN,y0,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, y0 is the initial condition of the
% ODE, and h is the stepsize.  You may also want to pass the function into
% the ODE the way |ode45| does (check lab 2).
%
% Note: you will need to use a loop to do this exercise.  
% You will also need to recall the Heun/Improved Euler algorithm learned in lectures.  
%

% Cleaning up workspace
close all; clear; clc;

%% Exercise 2
%
% Objective: Compare Heun with |ode45|.
%
% Specifics:  For the following initial-value problems (from lab 2, 
% exercises 1, 4-6), approximate the solutions with your function from
% exercise 1 (Improved Euler Method).
% Plot the graphs of your Improved Euler Approximation with the |ode45| 
% approximation.
%
% (a) |y' = y tan t + sin t, y(0) = -1/2| from |t = 0| to |t = pi|
%
% (b) |y' = 1 / y^2 , y(1) = 1| from |t=1| to |t=10|
%
% (c) |y' =  1 - t y / 2, y(0) = -1| from |t=0| to |t=10|
%
% (d) |y' = y^3 - t^2, y(0) = 1| from |t=0| to |t=1|
%
% Comment on any major differences, or the lack thereof. You do not need
% to reproduce all the code here. Simply make note of any differences for
% each of the four IVPs.

% For (a), the Improved Euler Method (IEM) and ode45 nearly identical
% solutions, except around pi/2. The exact solution is y(t) = -cos(t) / 2,
% which does not have any asymptotes. However, since the differential
% equation has tan(t) in its equation, when the IEM tries to calculate the
% slope near pi/2, it approaches negative or positive infinity. This
% behaviour is not representative of the actual analytical solution. On the
% other hand, ode45 appears to handle the discontinuity much better, as its
% solution appears smooth near pi/2.

% For (b), both the IEM and ode45 have very similar-looking plots. At a
% glance, there is no obvious difference between the two solutions.

% For (c), for the most part, the IEM and ode45 plots look nearly
% identical. However, when examining their behaviour around t = [1, 2], the
% IEM appears to have a smoother curve, whereas ode45 has some straight
% lines.

% For (d), neither the IEM nor ode45 have reasonable-looking solutions. The
% IEM solution shoots off to infinity around 0.8 and ode45 shoots off to
% infinity around 0.5.

% Cleaning up workspace
close all; clear; clc;

% (a)
% (a) |y' = y tan t + sin t, y(0) = -1/2| from |t = 0| to |t = pi|
f_a = @(t, y) y * tan(t) + sin(t);
t_0_a = 0;
t_N_a = pi;
y_0_a = -1/2;
h = 0.1;
[t, y] = improvedEuler(f_a, t_0_a, t_N_a, y_0_a, h);
plot(t, y, 'x-');
hold on;
soln_a = ode45(f_a, [t_0_a, t_N_a], y_0_a);
plot(soln_a.x, soln_a.y, 'x--');
title('a');
xlabel('t');
ylabel('y');
legend('Improved Euler Method', 'ode45', 'Location', 'Best');
snapnow;
clf;

% (b)
% (b) |y' = 1 / y^2 , y(1) = 1| from |t=1| to |t=10|
f_b = @(t, y) 1 / (y^2);
t_0_b = 1;
t_N_b = 10;
y_0_b = 1;
h = 0.1;
[t, y] = improvedEuler(f_b, t_0_b, t_N_b, y_0_b, h);
plot(t, y, 'x-');
hold on;
soln_b = ode45(f_b, [t_0_b, t_N_b], y_0_b);
plot(soln_b.x, soln_b.y, 'x--');
title('b');
xlabel('t');
ylabel('y');
legend('Improved Euler Method', 'ode45', 'Location', 'Best');
snapnow;
clf;

% (c)
% (c) |y' =  1 - t y / 2, y(0) = -1| from |t=0| to |t=10|
f_c  = @(t, y) 1 - (t * y) / 2;
t_0_c = 0;
t_N_c = 10;
y_0_c = -1;
h = 0.1;
[t, y] = improvedEuler(f_c, t_0_c, t_N_c, y_0_c, h);
plot(t, y, 'x-');
hold on;
soln_c = ode45(f_c, [t_0_c, t_N_c], y_0_c);
plot(soln_c.x, soln_c.y, 'x--');
title('c');
xlabel('t');
ylabel('y');
legend('Improved Euler Method', 'ode45', 'Location', 'Best');
snapnow;
clf;

% (d)
% (d) |y' = y^3 - t^2, y(0) = 1| from |t=0| to |t=1|
f_d = @(t, y) y^3 - t^2;
t_0_d = 0;
t_N_d = 1;
y_0_d = 1;
h = 0.1;
[t, y] = improvedEuler(f_d, t_0_d, t_N_d, y_0_d, h);
plot(t, y, 'x-');
hold on;
soln_d = ode45(f_d, [t_0_d, t_N_d], y_0_d);
plot(soln_d.x, soln_d.y, 'x--');
title('d');
xlabel('t');
ylabel('y');
legend('Improved Euler Method', 'ode45', 'Location', 'Best');
snapnow;
clf;

%% Exercise 3
%
% Objective: Use Euler's method and verify an estimate for the global error.
%
% Details: 
%
% (a) Use Euler's method (you can use
% euler.m from iode) to solve the IVP
%
% |y' = 2 t sqrt( 1 - y^2 )  ,  y(0) = 0|
%
% from |t=0| to |t=0.5|.
%
% (b) Calculate the solution of the IVP and evaluate it at |t=0.5|.
%
% (c) Read the attached derivation of an estimate of the global error for 
%     Euler's method. Type out the resulting bound for En here in
%     a comment. Define each variable.
%
% (d) Compute the error estimate for |t=0.5| and compare with the actual
% error.
%
% (e) Change the time step and compare the new error estimate with the
% actual error. Comment on how it confirms the order of Euler's method.

% Cleaning up workspace
close all; clear; clc;

% (a)
% Defining function and setting parameters
f = @(t, y) 2 * t * sqrt(1 - y^2);
t_0 = 0;
t_N = 0.5;
y_0 = 0;

% Solving IVP
t = linspace(t_0, t_N, 25);
y = euler(f, y_0, t);

% Plotting solution
plot(t, y, 'x--');
title('Solution to dy/dt = 2t * sqrt(1 - y^2)');
xlabel('t');
ylabel('y');
hold on;
f_exact = @(t) sin(t.^2);
[x, y] = plotExact(f_exact, t_0, t_N);
plot(x, y);
legend('Euler', 'Exact');

% (b)
% General solution: y = sin(t^2) + sin(C)
% Particular solution: y(t) = sin(t^2)
t = 0.5;
y_particular = sin(t^2);
fprintf('y(%g) = %g\n', t, y_particular);

% (c)
% E_n <= (1 + M)*dt / 2 * (e^(M*dt*n) - 1)
% E_n: error at step n
% M: upper bound of f, the partial derivative of f wrt y, and the partial
%    derivative of f wrt t between t = [0, 0.5]
% dt: step size
% n: step number

% (d)
M = 2;
dt = 0.01;

% Computing actual error
t = t_0:dt:t_N;
y = euler(f, y_0, t);
act_value = f_exact(0.5);
euler_value = y(length(y));
act_error = abs(act_value - euler_value);

% Computing error estimate
n = length(t) - 1;
err_func = @(M, dt, n) ((1 + M) * dt * (exp(M * dt * n) - 1)) / 2;
comp_error = err_func(M, dt, n);

fprintf('dt = %g\n', dt);
fprintf('Computed Error: %g\n', comp_error);
fprintf('Actual Error: %g\n', act_error);

% (e)
dt = 0.001;

% Computing actual error
t = t_0:dt:t_N;
y = euler(f, y_0, t);
f_exact = @(t) sin(t.^2);
act_value = f_exact(0.5);
euler_value = y(length(y));
act_error = abs(act_value - euler_value);

% Computing error estimate
n = length(t) - 1;
comp_error = err_func(M, dt, n);

fprintf('dt = %g\n', dt);
fprintf('Computed Error: %g\n', comp_error);
fprintf('Actual Error: %g\n', act_error);

% Confirms first order nature of Euler's method because decreasing the step
% size by a factor of 10 reduces the maximum possible error by a factor of
% 10 as well.

%% Adaptive Step Size
%
% As mentioned in lab 2, the step size in |ode45| is adapted to a
% specific error tolerance.
%
% The idea of adaptive step size is to change the step size |h| to a
% smaller number whenever the derivative of the solution changes quickly.
% This is done by evaluating f(t,y) and checking how it changes from one
% iteration to the next.

%% Exercise 4
%
% Objective: Create an Adaptive Euler method, with an adaptive step size |h|.
%
% Details: Create an m-file which accepts the variables |(t0,tN,y0,h)|, as 
% in exercise 1, where |h| is an initial step size. You may also want to 
% pass the function into the ODE the way |ode45| does.
%
% Create an implementation of Euler's method by modifying your solution to 
% exercise 1. Change it to include the following:
%
% (a) On each timestep, make two estimates of the value of the solution at
% the end of the timestep: |Y| from one Euler step of size |h| and |Z| 
% from two successive Euler steps of size |h/2|. The difference in these
% two values is an estimate for the error.
%
% (b) Let |tol=1e-8| and |D=Z-Y|. If |abs(D)<tol|, declare the step to be
% successful and set the new solution value to be |Z+D|. This value has
% local error |O(h^3)|. If |abs(D)>=tol|, reject this step and repeat it 
% with a new step size, from (c).
%
% (c) Update the step size as |h = 0.9*h*min(max(tol/abs(D),0.3),2)|.
%
% Comment on what the formula for updating the step size is attempting to
% achieve.

% The formula is trying to reduce the step size if it is determined that
% the value of D (i.e. the error) is too large. This ensures that the
% solution is within a given tolerance. The 0.9 is a safety factor to
% ensure success on the next step. The |min| and |max| functions prevent
% large changes from the previous stepsize.

%% Exercise 5
%
% Objective: Compare Euler to your Adaptive Euler method.
%
% Details: Consider the IVP from exercise 3.
%
% (a) Use Euler method to approximate the solution from |t=0| to |t=0.75|
% with |h=0.025|.
%
% (b) Use your Adaptive Euler method to approximate the solution from |t=0| 
% to |t=0.75| with initial |h=0.025|.
%
% (c) Plot both approximations together with the exact solution.

% Clearing workspace
close all; clear; clc;

% (a)
% Defining function and parameters
f = @(t,y) 2 * t * sqrt(1 - y^2);
t_0 = 0;
t_N = 0.75;
y_0 = 0;
h = 0.025;

t = t_0:h:t_N;
y = euler(f, y_0, t);
plot(t, y, 'x--');
hold on;

% (b)
[t, y] = adaptiveEuler(f, t_0, t_N, y_0, h);
plot(t, y);
hold on;

% (c)
f_exact = @(t) sin(t.^2);
[t, y] = plotExact(f_exact, t_0, t_N);
plot(t, y, '--');
title('Comparing Euler and Adaptive Euler Methods');
xlabel('t');
ylabel('y');
legend('Euler', 'Adaptive Euler', 'Exact', 'Location', 'Best');

%% Exercise 6
%
% Objective: Problems with Numerical Methods.
%
% Details: Consider the IVP from exercise 3 (and 5).
% 
% (a) From the two approximations calculated in exercise 5, which one is
% closer to the actual solution (done in 3.b)? Explain why.
% 
% (b) Plot the exact solution (from exercise 3.b), the Euler's 
% approximation (from exercise 3.a) and the adaptive Euler's approximation 
% (from exercise 5) from |t=0| to |t=1.5|.
%
% (c) Notice how the exact solution and the approximations become very
% different. Why is that? Write your answer as a comment.

% NOTE: Graph for Exercise 6 is at the bottom of the PDF.

% Clearing workspace
close all; clear; clc;

% (a)
% The Adaptive Euler Method is closer to the exact solution because the
% step size is constantly adjusted to fall within the given tolerance. On
% the other hand, the Euler Method only has a fixed step size, so if the
% function's derivative has a sudden change, then it will not be able to
% account for the rapid change as precisely.

% (b)
f = @(t,y) 2 * t * sqrt(1 - y^2);
t_0 = 0;
t_N = 1.5;
y_0 = 0;
h = 0.025;

f_exact = @(t) sin(t.^2);
[t, y] = plotExact(f_exact, t_0, t_N);
plot(t, y);
hold on;

t = t_0:h:t_N;
y = euler(f, y_0, t);
plot(t, y, 'x--');
hold on;

[t, y] = adaptiveEuler(f, t_0, t_N, y_0, h);
plot(t, y, '--');
title('Examining Various Numerical Methods');
xlabel('t');
ylabel('y');
legend('Exact', 'Euler', 'Adaptive Euler', 'Location', 'Best');
snapnow;

% (c)
% The exact solutions and approximations become very different around t = 1
% because the DE has sqrt(1 - y^2). When y > 1, the numerical methods try
% to compute the slope using the given DE, but get a result containing an
% imaginary component. As a result, the numerical methods become incapable
% of calculating solutions and the slopes stop changing.

%% General Functions
%
% Helper functions that I wrote and used throughout this lab.

function [x, y] = plotExact(f, a, b)
    % Allocating space for x, y vectors
    x = linspace(a, b);
    y = f(x);
end
