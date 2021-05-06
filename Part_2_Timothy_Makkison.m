% Where A is the stress tensor at a specific point in the coiled spring.
A = [0 0 0; 0 3.1 -1.4; 0 -1.4 4.5]
% B is a identity matrix
B = [1 0 0; 0 1 0; 0 0 1]

% Use function to calculate stress
max_stress = polyIteration(A,B,4,10,15)
% Print max principal stress
% msgbox(sprintf('The calculated maximum principal stress is %+2.1f MPa at the given point.',max_stress))



% Define the frequency distribution of the calculated maximum principal
% stresses.
princ_stress = 0.5:0.5:5
freq = [10, 15, 19, 24, 36, 50, 42, 29, 20, 5]
x1 = 0.5:0.1:5.5

% Plot frequency and stress data as a scatter chart.
plot(princ_stress,freq,'ko')

% lets additional graphs draw on top
hold on

natSpline = csape(princ_stress,freq,'second',[0,0])
y2 = ppval(natSpline,x1)
plot(x1,y2,'b-')

% Specify a range of values on the x and y axis
% Add padding to make values more readable.
xlim([0,5.5])
ylim([0,70])

% Assign title and labels.
title('Graph of maximum principal stress frequency.')
xlabel('Max. Principal Stress (MPa)')
ylabel('Stress frequency')

legend('Stress Frequency','Natural Spline')




% Calculate the area under the graph
area = trapezoidArea(freq,0.5)

% print area with 2 sf
% msgbox(sprintf('Area under the principal stress frequency graph is %+2.0f MPa calculated using the Trapezoidal rule.',area))

function areaSum = trapezoidArea(array,h)
    N = length(array);
    
%     initialize the value used to accumulate the area.
    areaSum = 0;
    
    %     Loop for every item in array -1 otherwise the final item wouldnt
    %     have an item to calulcate the trapezium with.
    for i = 1:1:N-1
%         Calculte the trapezoid using the trapezium rule and the next
%         value.
        area = (array(i) + array(i+1)) * (h/2);
%         Add area to total sum
        areaSum = areaSum + area;
    end
end

function lamb2 = polyIteration(A,B,lamb1,lamb2,iter)
% Uses polynomial iteration to calculate the max eigen vlaue value for A
% uses lamb1 and lamb2 as starting points
% iter defines the maximum number of runs before returning the best
% approximation.
    
    
%     Assert that starting lambda values are greater than 0
%     as eigen values cannot be below zero
    assert(lamb1 >= 0)
    assert(lamb2 >= 0)
    
    for i=1:1:iter
%       Break the for loop if values are the same aka have converged
%       This prevents the code dividing by zero and creating Nan
%       overwriting the previously correct value.
        if lamb1 == lamb2
            break
        end
%         Calculate Q
        Q1 = A - (B*lamb1);
        Q2 = A - (B * lamb2);
        
%         Calculate the determinant of each Q matrix
        detQ1 = det(Q1);
        detQ2 = det(Q2);
        
%         Calculate the next lambda estimate
        lamb3 = lamb2 - detQ2 * ((lamb1 - lamb2)/(detQ1 - detQ2));
        
%         Assign new estimate to lamb2 and assign lamb1 the old estimate
        lamb1 = lamb2;
        lamb2 = lamb3;
    end
end
