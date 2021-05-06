% Clear command window
clc
% Clear workspace to prevent external values altering the script.
clear

% Define suspension system equation

% Define specified parameters

% c is the damping coefficient in Ns/m
c = 1.5*10^3;
% k is the spring constant in N/m
k = 1.5*10^4;
% m is mass in kg
m = 3.6*10^3;

% Y0 is the initial displacement in metres
y0 = 0.3;

% Define calculated values
n = c/(2*m);
springConstPerMetre = k/m;
dampCoPerMass = c^2 / (4*m^2);

assert(springConstPerMetre > dampCoPerMass)

rootD = sqrt(springConstPerMetre - dampCoPerMass);

% Define symbol t for use in function
syms t
% Define the symbolic function of the car suspension.
systemFunc(t) = exp(-n*t) * y0 * (cos(rootD*t) + (n/rootD)*sin(rootD*t));




% Plot calculated data graph 

% define values for x
x1 = 0:0.01:10;

% Pass values of into function and plot
plot(x1,systemFunc(x1),'-b')

% Set x axis to strike through the origin
ax = gca;
ax.XAxisLocation ='origin'

% Set x and y labels
xlabel("Time (s)")
ylabel("Displacement (m)")

% Give graph a title
title("Graph of displacement in car suspension over time.")


% Calculate the first three roots

% Set minimum error value used to find root.
error = 1E-3

% Identify the first 3 roots as around 0.5, 2, 3.7
% Calculate roots iteratively using the Newton-Raphson second method.
root1 = double(findRoot(0.5, systemFunc, error))
root2 = double(findRoot(2, systemFunc, error))
root3 = double(findRoot(3.7, systemFunc, error))

% Print the first 3 roots.
msgbox(sprintf(['First three roots of the spring-damper system starting ' ...
    'at y0 = 0.3m are %2.2f, %2.2f and %2.2f seconds (s)' ...
    'calculated using the Newton-Raphson second method.'],root1,root2,root3))




% Comapare experimental and calculated data.

% Create the values for time.
time = 2.6:0.05:3.1

% Experimental displacement data
experimentDisp = [0.074, 0.094, 0.106, 0.113, 0.127, 0.145, 0.145, 0.148, 0.154, 0.162, 0.158]
% Expected displacement data calculated using the given function.
calculatedDisp = double(systemFunc(time))

[R,e] = calcErrorAndDetermination(experimentDisp,calculatedDisp)

msgbox(sprintf(['The values for coefficient of multiple determination and ' ...
    'standard error of the estimate were calculated to' ...
    ' be: \n R = %1.3f \n error = %1.3f'],R,e))





% Plot the experimental results data

figure

% Plot experimental displacemet data as a scatter graph.
plot(time,experimentDisp,'ko')

% Allow other plots to draw over the graph.
hold on

% Set x and y labels
xlabel("Time (s)")
ylabel("Displacement (m)")

% Give graph a title
title("Experimental data graph for the displacement of a car wheel over time.")

% Define an finer array of the same ranges as time.
t2 = 2.6:0.01:3.1;

% From theory saturation curve graph points are equal to 
% the inverse of a corresponding linear graph.
inverseTime = 1./time;
inverseExpDisp = 1./experimentDisp;

% Calculate the linear regression of the inverse values.
[a,b] = linear_regression(inverseTime,inverseExpDisp);

% Alpha is equal to the inverse of a.
Alpha = 1/a;
% Beta is equal to b/a.
Beta = b/a;

% Calculate the saturation curve using the equation (Alpha*x)/(Beta+x)
saturationCurve = (Alpha * t2)./(Beta+t2);

% Plot the saturation curve.
plot(t2,saturationCurve,'b')

% Adjust the range limits so all scatter points are visible.
xlim([2.5,3.2])
ylim([0.065,0.2])

% Define the legend of the graph.
legend('Displacement Data Points','Saturation Line Curve','Location',"southeast")


function [R,e] = calcErrorAndDetermination(y,g)
% Calculate the coefficient of multiple determination and 
% standard error of the estimate 
% Where y and g are arrays of the same length
% Where y is the found value and g is the predicted/calculated
% Returns arrays where R is the coefficient of multiple determination and
% e is the standard error of the estimate.

% Assert that inputs are the same size to prevent error and ensure
% a correct calculation.
    assert(length(y) == length(g))
    
    meanY = mean(y);
%     Calculate the total number of items
    N = length(y);
    
% Calculate the total sum of square and residual sum
    totalSumSquares = sum((g - meanY).^2)
    residualSumSquares =sum((y - meanY).^2)
    
% Calculate the coefficient of multiple determination
    R = (totalSumSquares/residualSumSquares)     
    
% Calculate the standard error of the estimate
    diffSum = sum((y-g).^2 )
    e = sqrt( diffSum/ (N-2) )
end    
function [a,b] = linear_regression(x,y)
% Calculate the linear rengression using two arrays.
    assert(length(x) == length(y))
    
% Calculate the sum of y,x,xy and x^2.
    x_sum = sum(x);
    y_sum = sum(y);
    xy_sum = sum(x.*y);
    x_square_sum = sum(x.^2);
    
% Get the length of the arrays
    N = length(x);
    
% To avoid rewriting code calculate the shared denominator.
    denom = N*x_square_sum - x_sum^2;
    
% Calculate a and b.
    a = (y_sum*x_square_sum - xy_sum*x_sum) / denom
    b = (N*xy_sum - y_sum*x_sum) / denom 

end
function x = findRoot(x0,func,error)
% Iteratively calculate the root to within a certain error using 
% the newton rhapson formula.

% Check required error is greater than 0
    assert(error > 0)
    
%     Find the first estimate for x
    x = newtonRhapson_2(x0,func);
%     Use the new value for x to find the height at that point
    y = func(x);
    
%     While height is not within the margin for error then loop
    while abs(y) > error
%         Calculate a new estimate for x
        x = double(newtonRhapson_2(x,func));
%         Calculate the height at that point
        y = double(func(x));
    end
end

function x = newtonRhapson_2(x0,func)
% Implementation of Newton Rhapson's second equation.
% Note it would be more efficicient to calculate the derivatives once and 
% inline this function into findRoot, however I believe this is easier to
% read.

% Calulate the first derivative of the given function
    derFunc = diff(func);
% Calculate the second derivative of the given function
    secondDerFunc  = diff(func,2);
    
% Calculate the slope when x = x0
    slope = derFunc(x0);
% Calculate the derivative of the slope when x = x0
    derSlope = secondDerFunc(x0);
    
    x = x0 + (derSlope/(2*slope) - slope/func(x0))^-1;
end
