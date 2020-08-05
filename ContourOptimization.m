clear all; clc; close all
% From Carsim paper
% At Rho = 30 m -> U = 1.91
% At Rho = 60 m -> U = 1.95
y1 = linspace(-10,10); %Resonable Angles
y2 = linspace(0,40); % Reasonable Speed Ranges (m/s)
% 40 m/s ~ 90 mph ~ 144 km/hr
L = 2.5; % m
g = 9.81; % m/s^2
%AASHTO values
e = 6; mu = .4;
U  = 1.95; K = 1/60;
[X1,X2] = meshgrid(y1,y2);
Z = X1 - (53.7*L+U*X2.^2)*K;
figure; meshc(X1,X2,Z); hold on
Z2 = - X2.^2*K/g + (mu + 0.01*e)/(1-0.01*mu*e);
meshc(X1,X2,Z2); colorbar; 
xlabel('y_1'); ylabel('y_2')
