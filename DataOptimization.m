clear; close all; clc
load('CVF9LatX.mat')
load('CVF9LongY.mat')
x2 = LatX'; y2 = LongY';   
%x2(1:2:end) = [];  y2(1:2:end) = [];
x2 = unique(x2);
y2 = unique(y2);
x2 = x2(1:numel(y2));X = [x2',y2'];
[L2,R2,K2] = curvature(X);
figure(1); plot(L2,R2); grid on;
title('Curvature radius \rho vs. Cumulative curve length')
%The Radius of Curvature is High at the Ends and Small at the middle
%Which is opposite for the Curvature Kappa
xlabel('Length of Road'); ylabel('Radius \rho')
figure(2);
h = plot(x2,y2); grid on; axis equal; set(h,'marker','.');
xlabel('X Coordinate'); ylabel('Y Coordinate')
title('Road with Curvature Vectors')
hold on
quiver(x2',y2',K2(:,1),K2(:,2)); hold off
% --- Not Used
% figure(3); hold on;
% plot(x2,sqrt(K2(:,2).^2+K2(:,1).^2))
% xlabel('X-Coordinate'); ylabel('Curvature \kappa'); grid on;
% title('Curvature \kappa vs X Coordinate')
% figure(4); hold on;
% plot(y2,sqrt(K2(:,2).^2+K2(:,1).^2))
% xlabel('Y-Coordinate'); ylabel('Curvature \kappa'); grid on;
% title('Curvature \kappa vs Y Coordinate')

KK = sqrt(K2(:,2).^2+K2(:,1).^2);
figure(5); scatter(L2,KK); grid on
xlabel('Lenght of Road'); ylabel('Curvature \kappa')
title('Curvature \kappa vs. Cumulative curve length')

%  Smoothing Technique on Curvature----------------
close all
figure(1050)
x = L2; y = KK;
yy1 = smooth(x,y,0.15,'loess');  %Span of 15%
yy2 = smooth(x,y,0.20,'rloess');

%subplot(2,1,1)
plot(x,y,'b.',x,yy1,'r-'); ylim([0,.3]); grid on
legend('Original data','Smoothed data using ''loess''',...
       'Location','best')
xlabel('Segment S (m)'); ylabel('Curvature K')
%subplot(2,1,2)
figure(1051)
plot(x,y,'b.',x,yy2,'r-'); grid on;  ylim([0,.3])
xlabel('Segment S (m)'); ylabel('Curvature K')
legend('Original data','Smoothed Data',...
       'Location','best')
   
% Adding Spline for Smoothed Function
figure(1055)

s = L2; k = KK; dt = 0.01;
n = numel(s);
[A1] = Spliny(s,k,n,dt);
% This is the curvature, without being smoothed, and 
% a spline interpolation has been performed to obtain all of its 
% constituent coefficients. 


% From data Figure 1051 rloess 20%
figure(1056)

s = x; k = yy2; dt = 0.01;
n = numel(s);
[A2] = Spliny(s,k,n,dt);
% This is the curvature, after being smoothed, and 
% a spline interpolation has been performed to obtain all of its 
% constituent coefficients. 

% NOTE: Matrix is close to singular/badly scaled.




