clc; clear all; close all
% Same as CurveFitlsq but uses a higher amount of data points
% Behavior is good because even though there should be 
% a "drop" the function compensates and does not generate a drop
% The only problem is that the function only behaves properly when the data
% starts at 0 for y-axis
% Also, Here I am trying different levels of gaussian noise
s = 1:.01:25;
%y = [0 2 4.8 5.2 5 5.6];
n = numel(s)-1;
y1 = 2.*s(1:n/2) - 1;
y2 = 23*ones(1,n/2);

for k = 1:9

y1o = awgn(y1,k,'measured');
y2o = awgn(y2,k,'measured');

y = [y1o y2o];
x0 = [0.9 2.9 4.89 7.85 11.50];

fun1 = @(x,s) ((x(5)./(x(2)-x(1))).*(s - x(1))).*(heaviside(s-x(1)) - heaviside(s-x(2))) +...
     x(5).*(heaviside(s-x(2))-heaviside(s-x(3))) + ...
( ( x(5)./(x(4)-x(3))).*(-s+x(3))+ x(5) ).*(heaviside(s-x(3)) - heaviside(s-x(4))); 

x = lsqcurvefit(fun1,x0,s(1:end-1),y)
times = linspace(s(1),s(end-1));
figure
hold on; plot(s(1:end-1),y,'bo')
plot(times,fun1(x,times),'k-','linewidth',2)
xlim([times(1), times(end)+5])
legend('Data','Fitted Response','location','best'); 
title('Data and Fitted Curve'); grid on

end
