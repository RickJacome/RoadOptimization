clear; close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this example, no noise is added because is the GPS trial
% (noise should be part of the data)
% Smoothing is performed to the Curvature Data coming from the
% noisy x,y coordinates, not to the data itself
% Recall that for filtering directions (vectors), 
%the smoothing did an excellent job (on GPSCVCXY.mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GPS DATA
load('CVF9LatX.mat'); load('CVF9LongY.mat');
%Ideal AASHTO
%load('IdealXm.mat'); load('IdealYm.mat');
x2 = LatX'; y2 = LongY';   
%x2 = xm'; y2 = ym';
x2 = unique(x2); y2 = unique(y2);
x2 = x2(1:numel(y2));
X = [x2',y2'];
[L,R,K] = curvature(X);
K(1,:) = []; K(end,:) = []; L(1,:) = []; L(end,:) = [];
x2(1) = []; x2(end) = []; y2(1) = []; y2(end) = [];
figure; plot(x2,y2);
xlabel('X Coordinate (m)'); ylabel('Y Coordinate (m)')
title('Raw Road Data')
figure;
h = plot(x2,y2); grid on; axis equal; set(h,'marker','.');
xlabel('X Coordinate (m)'); ylabel('Y Coordinate (m)')
title('Road with Curvature Vectors')
hold on
quiver(x2',y2',K(:,1),K(:,2)); hold off  
y = sqrt(K(:,1).^2 + K(:,2).^2);
s = L;
figure; plot(s,y)
yy1 = smooth(s,y,0.2,'loess');
s = unique(s); yy1 = unique(yy1);
s = s(1:numel(yy1));
y = yy1;

% Initial Conditions, NEVER repeat them.
x0 = [100 200 300 400 500];
% Curvature Model M.1
M1 = @(x,s) ((x(5)./(x(2)-x(1))).*(s - x(1))).*(heaviside(s-x(1)) - heaviside(s-x(2))) +...
     x(5).*(heaviside(s-x(2))-heaviside(s-x(3))) + ...
( ( x(5)./(x(4)-x(3))).*(-s+x(3))+ x(5) ).*(heaviside(s-x(3)) - heaviside(s-x(4))); 
% Pr.1
fprintf('Pr. 1, Least Squares Min. Has finalized');
options = optimset('Display','off');
x = lsqcurvefit(M1,x0,s(1:end),y,[],[],options)
snew = linspace(s(1),s(end),100); % <--- This defines the 
% size of the "K_vector".
figure; hold on; 
plot(s,y,'bo');
xlabel('S-Segment (m)'); ylabel ('Curvature(m^{-1})');
plot(snew,M1(x,snew),'k-','linewidth',2);
xlim([snew(1), snew(end)+5]);
legend('Data','Fitted Response','location','best'); 
title('Data and Fitted Curve'); grid on

% -------------------------
%Parameters 
global K_temp e g mu U
% Vehicle Only
L = 2.5;  %U = 1.95;
U = 3;
% Road Only
e = 6; mu = 0.3;
% Both
g = 9.81; K_vector = M1(x,snew);
% -------------------------
%Iterative Optimization Routine for Pr.2 given Optimized M.1
for i = 1:length(K_vector)
K_temp = K_vector(i);  
% Objective Function Pr.2
fun = @(x)  x(1) - (53.7*L + U*x(2)^2/g)*K_temp;    
%C.1 (Bounds)
% lb = [-3,25]; 
% ub = [3,60]; 
lb = [-3,25]; % -3 < x1 < 3;
ub = [3,35];  % 55 < x2 < 80; mph
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
fprintf('Pr. 2 Has finalized \n');
vOpt = Op(:,2);
figure; plot(snew,vOpt,'linewidth',3)
ylim([max(vOpt)-10 max(vOpt)+10])
title('Segment Length vs Velocity Optimized'); grid on
xlabel('S-Segment (m)'); ylabel ('Velocity (m/s)');
figure; plot(M1(x,snew),vOpt)
ylim([max(vOpt)-10 max(vOpt)+10])
title('Curvature vs Velocity Optimized'); grid on;
xlabel('S-Segment (m)'); ylabel ('Curvature(m^{-1})');

% Nonlinear Constaints (Not bounds)
function [c,ceq] = EqConstraint(x)
global K_temp e g mu 
%Pr.2
% Nonlinear Inequality Constraints
c = x(2)^2*K_temp/g - (mu + 0.01*e)/(1-0.01*mu*e);
% Nonlinear Equality Constraints
ceq = [];
end


