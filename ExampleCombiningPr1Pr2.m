clc; clear all; close all
% This code takes a curvature model (M.X) and
% implements both Least Squares Optimization and  The Nonlinear Iterative
% Optimization together.
% It creates a loop of optimization problems related to the Curvature
% model provided with Pr. 2. Such that an optimal velocity is found for 
% every curvature datapoint for any general [s] road. 
%------------------------------
% Road Curvature Decomposition Data
s = 1:.01:25; n = numel(s)-1;
y1 = (2.*s(1:n/2) - 3)*1e-3;
y2 = 26*ones(1,n/2)*1e-3;
% Added Gaussian noise.
y1o = awgn(y1,25,'measured'); y2o = awgn(y2,25,'measured');
y = [y1o y2o];
% Initial Conditions, NEVER repeat them.
x0 = [0.9 2.9 4.89 7.85 11.50];
% Curvature Model M.1
M1 = @(x,s) ((x(5)./(x(2)-x(1))).*(s - x(1))).*(heaviside(s-x(1)) - heaviside(s-x(2))) +...
     x(5).*(heaviside(s-x(2))-heaviside(s-x(3))) + ...
( ( x(5)./(x(4)-x(3))).*(-s+x(3))+ x(5) ).*(heaviside(s-x(3)) - heaviside(s-x(4))); 
% Pr.1
fprintf('Pr. 1, Least Squares Min. Has finalized \n');
options = optimset('Display','off');
x = lsqcurvefit(M1,x0,s(1:end-1),y,[],[],options)
snew = linspace(s(1),s(end-1),100); % <--- This defines the 
% size of the "K_vector".
figure; hold on; 
plot(s(1:end-1),y,'bo');
plot(snew,M1(x,snew),'k-','linewidth',2);
xlim([snew(1), snew(end)+5]);
legend('Data','Fitted Response','location','best'); 
title('Data and Fitted Curve'); grid on

% -------------------------
%Parameters 
global K_temp e g mu U
% Vehicle Only
L = 2.5;  U = 1.95;
% Road Only
e = 6; mu = 0.4;
% Both
g = 9.81; K_vector = M1(x,snew);
% -------------------------
%Iterative Optimization Routine for Pr.2 given Optimized M.1
for i = 1:length(K_vector)
K_temp = K_vector(i);  
% Objective Function Pr.2
fun = @(x)  x(1) - (53.7*L + U*x(2)^2/g)*K_temp;    
%C.1 (Bounds)
lb = [-3,10]; % -3 < x1 < 3;
ub = [3,50];  % 0 < x2 < 60;
% There are no linear constraints, so set those arguments to |[]|. 
A = [];  b = []; % Linear In-equality Constraints
Aeq = []; beq = [];  % Linear Equality Constraints
%Initial Conditions
x0 = [1/4,1/4];  
%Constraints as an annoynomous function
nonlcon = @EqConstraint;
options = optimoptions('fmincon','Display','off');
Op(i,:) = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);  
end
fprintf('Pr. 2 Nonlinear Optimization Has finalized \n');
figure; plot(snew,Op(:,2))
title('Segment Length vs Velocity Optimized'); grid on
figure; plot(M1(x,snew),Op(:,2))
title('Curvature vs Velocity Optimized'); grid on

% Nonlinear Constaints (Not bounds)
function [c,ceq] = EqConstraint(x)
global K_temp e g mu 
%Pr.2
% Nonlinear Inequality Constraints
c = x(2)^2*K_temp/g - (mu + 0.01*e)/(1-0.01*mu*e);
% Nonlinear Equality Constraints
ceq = [];
end

