clc; clear all; close all
s = 1:.01:12; n = numel(s)-1;
%y = [0 2 4.8 5.2 5 5.6];
y1 = 2.*s(1:n/2) - 1;   y2 = 11*ones(1,n/2);
y1o = awgn(y1,20,'measured'); y2o = awgn(y2,20,'measured');

y = [y1o y2o];
x0 = [0.9 2.9 6 max(y)];
% Find the "minimized error".
% fun1 = @(x,s) ((x(4)./(x(2)-x(1))).*s - x(1).*x(4)./(x(2)-x(1))).*(heaviside(s-x(1)) - heaviside(s-x(2))) +...
%      x(4).*(heaviside(s-x(2))-heaviside(s-x(3)));
fun1 = @(x,s) ((x(4)./(x(2)-x(1))).*(s - x(1))).*(heaviside(s-x(1)) - heaviside(s-x(2))) +...
     x(4).*(heaviside(s-x(2))-heaviside(s-x(3))); 
x = lsqcurvefit(fun1,x0,s(1:end-1),y)
times = linspace(s(1),s(end-1));
hold on; plot(s(1:end-1),y,'bo')
plot(times,fun1(x,times),'k-','linewidth',2)
legend('Data','Fitted Response'); title('Data and Fitted Curve')
xlim([times(1), times(end)+5])

x0 = [0.9 2.9 6 10 max(y)];
fun1 = @(x,s) ((x(5)./(x(2)-x(1))).*(s - x(1))).*(heaviside(s-x(1)) - heaviside(s-x(2))) +...
     x(5).*(heaviside(s-x(2))-heaviside(s-x(3))) + ...
( ( x(5)./(x(4)-x(3))).*(-s+x(3))+ x(5) ).*(heaviside(s-x(3)) - heaviside(s-x(4))); 

x = lsqcurvefit(fun1,x0,s(1:end-1),y)
times = linspace(s(1),s(end-1));
figure
hold on; plot(s(1:end-1),y,'bo')
plot(times,fun1(x,times),'k-','linewidth',2)
legend('Data','Fitted Response'); title('Data and Fitted Curve')
xlim([times(1), times(end)+5])

