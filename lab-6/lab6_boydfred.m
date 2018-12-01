%% Laplace Transform Lab: Solving ODEs using Laplace Transform in MATLAB
%
% This lab will teach you to solve ODEs using a built in MATLAB Laplace 
% transform function |laplace|. Also in this lab, you will write your own
% ODE solver using Laplace transforms and check whether the result yields
% the correct answer.
%
% You will learn how to use the |laplace| routine. 
% 
% There are five (5) exercises in this lab that are to be handed in.  
% Write your solutions in the template, including appropriate descriptions 
% in each step. Save the m-file and submit it on Quercus.
%
% Include your name and student number in the submitted file.
%
% MAT292, Fall 2018, Stinchcombe & Khovanskii, modified from
% MAT292, Fall 2017, Stinchcombe & Sinnamon, modified from
% MAT292, Fall 2015, Sousa, based on 
% MAT292, Fall 2013, Sinnamon & Sousa

%% Student Information
%
%  Student Name: Frederick Boyd
%
%  Student Number:
%

%% Using symbolic variables to define functions
% 
% Recall the use of symbolic variables and function explained in the MATLAB
% assignment #2.

syms t s x y

f = cos(t)
h = exp(2*x)


%% Laplace transform and its inverse

% The routine |laplace| computes the Laplace transform of a function

F=laplace(f)

%%
% By default it uses the variable |s| for the Laplace transform
% But we can specify which variable we want:

H=laplace(h)
laplace(h,y)

% Observe that the results are identical: one in the variable |s| and the
% other in the variable |y|

%% 
% We can also specify which variable to use to compute the Laplace
% transform:

j = exp(x*t)
laplace(j)
laplace(j,x,s)

% By default, MATLAB assumes that the Laplace transform is to be computed
% using the variable |t|, unless we specify that we should use the variable
% |x|

%% 
% We can also use inline functions with |laplace|. When using inline
% functions, we always have to specify the variable of the function.

l = @(t) t^2+t+1
laplace(l(t))

%% 
% MATLAB also has the routine |ilaplace| to compute the inverse Laplace
% transform

ilaplace(F)
ilaplace(H)
ilaplace(laplace(f))

%% 
% If |laplace| cannot compute the Laplace transform, it returns an
% unevaluated call.

g = 1/sqrt(t^2+1)
G = laplace(g)

%% 
% But MATLAB "knows" that it is supposed to be a Laplace transform of a
% function. So if we compute the inverse Laplace transform, we obtain the
% original function

ilaplace(G)

%%
% The Laplace transform of a function is related to the Laplace transform 
% of its derivative:

syms g(t)
laplace(diff(g,t),t,s)


%% Exercise 1
%
% Objective: Compute the Laplace transform and use it to show that MATLAB
% 'knows' some of its properties.
%
% Details:  
%
% (a) Define the function |f(t)=exp(2t)*t^3|, and compute its Laplace
%   transform |F(s)|.
% (b) Find a function |f(t)| such that its Laplace transform is
%   |(s - 1)*(s - 2))/(s*(s + 2)*(s - 3)|
% (c) Show that MATLAB 'knows' that if |F(s)| is the Laplace transform of
%   |f(t)|, then the Laplace transform of |exp(at)f(t)| is |F(s-a)|
% 
% (in your answer, explain part (c) using comments).      
%
% Observe that MATLAB splits the rational function automatically when
% solving the inverse Laplace transform.

% Clearing workspace...
close all; clear; clc;

% Defining symbolic variables
syms x s t

% (a)
f = @(t) exp(2*t) * t^3;
F = laplace(f(t));
disp(F);

% (b)
G = @(s) ((s - 1)*(s - 2))/(s*(s + 2)*(s - 3));
g = ilaplace(G(s));
disp(g);

% (c)
syms a f(t)
F_1 = laplace(f(t));
F_2 = laplace(exp(a*t) * f(t));
disp(F_1);
disp(F_2);

% MATLAB knows that the laplace transform of exp(at)*f(t) is F(s-a) because
% when calculating the laplace transform of exp(at)*f(t), the result
% appears very similar. The only difference between the two laplace
% transforms is a shift of F(s-a) when f(t) is multiplied by exp(at).

%% Heaviside and Dirac functions
%
% These two functions are builtin to MATLAB: |heaviside| is the Heaviside
% function |u_0(t)| at |0|
%
% To define |u_2(t)|, we need to write

f=heaviside(t-2)
ezplot(f,[-1,5])

% The Dirac delta function (at |0|) is also defined with the routine |dirac|

g = dirac(t-3)

% MATLAB "knows" how to compute the Laplace transform of these functions

laplace(f)
laplace(g)


%% Exercise 2
%
% Objective: Find a formula comparing the Laplace transform of a 
%   translation of |f(t)| by |t-a| with the Laplace transform of |f(t)|
%
% Details:  
%
% * Give a value to |a|
% * Let |G(s)| be the Laplace transform of |g(t)=u_a(t)f(t-a)| 
%   and |F(s)| is the Laplace transform of |f(t)|, then find a 
%   formula relating |G(s)| and |F(s)|
%
% In your answer, explain the 'proof' using comments.

% Clearing workspace...
close all; clear; clc;

% Defining symbolic variables
syms f(t) s t a

% Defining constant functions
% f = @(t) exp(2*t) * t^3;

% Calculating laplace transforms using a for loop with various values of
% |a|
trials = 3;
for i = 1:trials
    % Assigning values and defining functions
    a = i;
    g = @(t) heaviside(t-a)*f(t-a);
    F = laplace(f(t));
    G = laplace(g(t));
    % Displaying results
    fprintf("a = %g\n", i);
    disp(F);
    disp(G);
end

% G(s) = F(s) * exp(-a*s)
% This "proof" can be seen by performing multiple Laplace transformations.
% The only difference between F(s) and G(s) is the multiplication by an
% exponential function in G(s). As the value of |a| increases, the constant
% in the exponential function increases as well.

%% Solving IVPs using Laplace transforms
%
% Consider the following IVP, |y''-3y = 5t| with the initial
% conditions |y(0)=1| and |y'(0)=2|.
% We can use MATLAB to solve this problem using Laplace transforms:

% First we define the unknown function and its variable and the Laplace
% tranform of the unknown

syms y(t) t Y s

% Then we define the ODE

ODE=diff(y(t),t,2)-3*y(t)-5*t == 0

% Now we compute the Laplace transform of the ODE.

L_ODE = laplace(ODE)

% Use the initial conditions

L_ODE=subs(L_ODE,y(0),1)
L_ODE=subs(L_ODE,subs(diff(y(t), t), t, 0),2)

% We then need to factor out the Laplace transform of |y(t)|

L_ODE = subs(L_ODE,laplace(y(t), t, s), Y)
Y=solve(L_ODE,Y)

% We now need to use the inverse Laplace transform to obtain the solution
% to the original IVP

y = ilaplace(Y)

% We can plot the solution

ezplot(y,[0,20])

% We can check that this is indeed the solution

diff(y,t,2)-3*y


%% Exercise 3
%
% Objective: Solve an IVP using the Laplace transform
%
% Details: Explain your steps using comments
%
%
% * Solve the IVP
% *   |y'''+2y''+y'+2*y=-cos(t)
% *   |y(0)=0|, |y'(0)=0|, and |y''(0)=0|
% * for |t| in |[0,10*pi]|
% * Is there an initial condition for which |y| remains bounded as |t| goes to infinity? If so, find it.

% Clearing workspace
close all; clear; clc;

% Configuring symbolic variables
syms y(t) t Y s

% Defining ODE and initial conditions
y_diff_3 = diff(y(t), t, 3);
y_diff_2 = diff(y(t), t, 2);
y_diff_1 = diff(y(t), t, 1);
ODE = y_diff_3 + 2 * y_diff_2 + y_diff_1 + 2 * y(t) == -cos(t);
y_0_0 = 0;
y_0_1 = 0;
y_0_2 = 0;

% Calculating Laplace transform and substituting initial conditions
L_ODE = laplace(ODE);
L_ODE = subs(L_ODE, y(0), y_0_0);
L_ODE = subs(L_ODE, subs(y_diff_1, t, 0), y_0_1);
L_ODE = subs(L_ODE, subs(y_diff_2, t, 0), y_0_2);

% Factoring out Y
L_ODE = subs(L_ODE, laplace(y(t), t, s), Y);
Y = solve(L_ODE, Y);

% Solving ODE using the inverse Laplace transform
y = ilaplace(Y);

% Plotting function
ezplot(y, [0, 10*pi]);

% There is no initial condition for which |y| remains bounded as |t| goes
% to infinity.

%% Exercise 4
%
% Objective: Solve an IVP using the Laplace transform
%
% Details:  
% 
% * Define 
% *   |g(t) = 3 if 0 < t < 2|
% *   |g(t) = t+1 if 2 < t < 5|
% *   |g(t) = 5 if t > 5|
%
% * Solve the IVP
% *   |y''+2y'+5y=g(t)|
% *   |y(0)=2 and y'(0)=1|
%
% * Plot the solution for |t| in |[0,12]| and |y| in |[0,2.25]|.
%
% In your answer, explain your steps using comments.

% Clearing workspace
close all; clear; clc;

% Configuring symbolic variables
syms y(t) t Y s

% Defining functions and initial conditions
% g(t) was calculated and simplified by hand
g = @(t) 3*heaviside(t) + (t-2) * heaviside(t-2) + (-t+4) * heaviside(t-5);
y_diff_2 = diff(y(t), t, 2);
y_diff_1 = diff(y(t), t, 1);
ODE = y_diff_2 + 2 * y_diff_1 + 5 * y(t) == g(t);
y_0_0 = 2;
y_0_1 = 1;

% Calculating Laplace transform and substituting inital conditions
L_ODE = laplace(ODE);
L_ODE = subs(L_ODE, y(0), y_0_0);
L_ODE = subs(L_ODE, subs(y_diff_1, t, 0), y_0_1);

% Factoring out Y
L_ODE = subs(L_ODE, laplace(y(t), t, s), Y);
Y = solve(L_ODE, Y);

% Solving ODE using the inverse Laplace transform
y = ilaplace(Y);

% Plotting function
t = linspace(0, 12, 250);
y_plot = subs(y);
plot(t, y_plot);
ylim([0, 2.25]);
title('Solution to Piecewise ODE');
xlabel('t');
ylabel('y');

%% Exercise 5
%
% Objective: Solve an IVP with a periodic function using the Laplace
% transform.
%
% Details:
%
% * Follow the instructions on the website
%  http://instruct.math.lsa.umich.edu/lecturedemos/ma216/docs/7_5/
%
% * Do the part labelled |Outside of Lecture|
%
% * Note that |u(t-a)| is the Heaviside function |u_a(t)| defined 
% in our textbook
%
% * Check Theorem 5.5.3 (page 349) to know how to define the Laplace 
% transform of a periodic function like the one in this exercise (and check
% the function |int| on MATLAB for symbolic integration).
%
% * Hint: Use only the second-order DE given on the linked website
% * Hint: You should use the numbers provided (the less precise ones)
% * Hint: You should obtain an exact solution. In the process, you should need symbolic integration(s).

% Clearing workspace
close all; clear; clc;

% Configuring symbolic variables
syms y(t) t Y s i;

% Defining functions, constants, and initial conditions
k1 = 0.12;
k2 = 0.38;
k3 = 0.04;
I0 = 24;   % Magnitude of dosing
t0 = 6;    % Length of dosing
t1 = 24;   % Time between doses
y_0_0 = 0;
y_0_1 = 0;

% Finding Laplace transform of the ODE's left side
y_diff_2 = diff(y(t), t, 2);
y_diff_1 = diff(y(t), t, 1);
ODE_LEFT = y_diff_2 + (k1 + k2 + k3) * y_diff_1 + k3 * (k1 + k2) * y(t);
L_ODE_LEFT = laplace(ODE_LEFT);

% Finding Laplace transform of the ODE's right side
% Using theorem 5.5.3 from the textbook to find Laplace transform of
% periodic functions
f(t) = (I0 / t0) * (heaviside(t) - heaviside(t - t0));
F(t) = int(exp(-s * t) * f(t), t, [0, t1]) / (1 - exp(-s * t1));
% Exploit linear properties of the Laplace transform by factoring out k1
L_ODE_RIGHT = k1 * F(t);
L_ODE = L_ODE_LEFT == L_ODE_RIGHT;

% Substituting initial conditions into ODE
L_ODE = subs(L_ODE, y(0), y_0_0);
L_ODE = subs(L_ODE, subs(y_diff_1, t, 0), y_0_1);

% Factoring out Y
L_ODE = subs(L_ODE, laplace(y(t), t, s), Y);
Y = solve(L_ODE, Y);

% Solving ODE using the inverse Laplace transform
y = ilaplace(Y);
fprintf("Exact Solution: %s\n", y);

% Plotting solution
t = linspace(0, 250, 500);
y_plot = subs(y);
plot(t, y_plot);
title('Antihistamine Model with Repeated Doses');
xlabel('t');
ylabel('y');
