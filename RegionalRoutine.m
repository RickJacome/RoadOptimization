%%%%---------------------------------------------
% 2-28-2021
% For this SemiDynamic Routine to work
% (1) Select your Data i.e. GPS, Google Earth, AASHTO, etc.
% (2) Select a range that occupies one "curve segment" i.e.
% when does curvature change considerably. 
% (3) Select appropiate Initial Conditions i.e. proportional to 
% segment length.
% External files used: curvature.m direction.m
%%%---------------------------------------------
clear; close all; clc
%%%---------------------------------------------
%(1)
%Google Earth Data
% load('GPS1Xft.mat'); load('GPS1Yft.mat');
% x2 = GPSX; y2 = GPSY;
% x2 = x2'*.3048; y2 = y2'*.3048;
% % Normalization (if needed)-----
% x2 = x2-min(x2);
% y2 = y2-min(y2);
% x2 = x2(100:180);
% y2 = y2(100:180);
% --------------------------
%GPS DATA
% load('CVF9LatX.mat'); load('CVF9LongY.mat');
% x2 = LatX'; y2 = LongY'; 
%Ideal AASHTO
load('MichXm.mat'); load('MichYm.mat');  
x2 = xm'; y2 = ym';
%%%---------------------------------------------
x2 = unique(x2,'stable'); y2 = unique(y2,'stable');
x2 = x2(1:numel(y2));
X = [x2',y2'];
[L,R,K] = curvature(X);
K(1,:) = []; K(end,:) = [];
L(1,:) = []; L(end,:) = [];
xlabel('Length of Road (m)'); ylabel('Radius \rho (m)')
figure(1);
x2(1) = []; x2(end) = [];
y2(1) = []; y2(end) = [];
h = plot(x2,y2); grid on; axis equal; set(h,'marker','.');
xlabel('X Coordinate (m)'); ylabel('Y Coordinate (m)')
title('Road with Curvature Vectors')
hold on
quiver(x2',y2',K(:,1),K(:,2)); hold off  
% ------------
y = sqrt(K(:,1).^2 + K(:,2).^2);
s = L;
%-------
[O1,O2] = direction(K);
e1 = cosd(O2); e2 = sind(O2);
figure(200);hold on; h1 = plot(x2,y2); grid on; axis equal; set(h1,'marker','.','Linewidth',3);
quiver(x2',y2',e1,e2); hold off
%title('Road with Velocity Vectors')
xlabel('X Coordinate (m)'); ylabel('Y Coordinate (m)');
figure; plot(s,y)
%%%---------------------------------------------
%(2)
%%% Road Section-----
ni = 1;
ne = numel(x2);
%%%---------------------------------------------
figure; plot(x2(ni:ne),y2(ni:ne),'linewidth',2); grid on;
xlabel('X Coordinate (m)'); ylabel('Y Coordinate (m)');
title('A')
figure; plot(s(ni:ne),y(ni:ne));
grid on;xlabel('Segment Length s (m)'); ylabel('Curvature \kappa (m^{-1})')
title('B')
ySmoo = smooth(s(ni:ne),y(ni:ne),0.15,'loess');
%ySmoo = y(ni:ne);
ySmoo = ySmoo(ni:10:ne);
sSmoo = s(ni:10:ne);
figure; plot(sSmoo,ySmoo); grid on;
xlabel('Segment Length s (m)'); ylabel('Curvature \kappa (m^{-1})')
%---------------------------------------------------------------
%(3)
%From this point forward, I am doing the Optimization Semi-Dynamic Routine
% Initial Conditions, NEVER repeat them.
%Ideal AASHTO IC.
%x0 = [ 100 200 300 400 max(ySmoo)];
x0 = [sSmoo(1) mean(sSmoo) 1.10*mean(sSmoo) sSmoo(end) max(ySmoo)];
%Google Earth IC.
%Constant Velocity
%x0 = [sSmoo(1) .70*mean(sSmoo) 1.25*mean(sSmoo) .90*sSmoo(end) max(ySmoo)];
% Recommended values for non-normalized data
% x0 = [3200 3500 4000 4200 max(ySmoo)];
% Curvature Model M.1
M1 = @(x,s) ((x(5)./(x(2)-x(1))).*(s - x(1))).*(heaviside(s-x(1)) - heaviside(s-x(2))) +...
     x(5).*(heaviside(s-x(2))-heaviside(s-x(3))) + ...
( ( x(5)./(x(4)-x(3))).*(-s+x(3))+ x(5) ).*(heaviside(s-x(3)) - heaviside(s-x(4))); 
% Pr.1
fprintf('Pr. 1, Least Squares Min. Has finalized');
options = optimset('Display','off');


x = lsqcurvefit(M1,x0,sSmoo,ySmoo,[],[],options)
snew = linspace(sSmoo(1),sSmoo(end),100); % <--- This defines the 
% size of the "K_vector".
figure; hold on; 
plot(sSmoo,ySmoo,'bo');
plot(snew,M1(x,snew),'k-','linewidth',2);
xlim([snew(1), snew(end)+5]);
legend('Data','Fitted Response','location','best'); 
title('Data and Fitted Curve'); grid on
xlabel('Segment Length s (m)'); ylabel('Curvature \kappa (m^{-1})')

% -------------------------
%Parameters 
global K_temp e g mu
% Road Only
%e = 12; mu = 0.4;
%e = 4; mu = 0.9;
e = 4; mu = 0.1;
%e = 4; mu = 0.4;
%e = 4; mu = 0.5;
%e = 4; mu = 0.8;
% Both
g = 9.81; K_vector = M1(x,snew);
% -------------------------
%Iterative Optimization Routine for Pr.2 given Optimized M.1
for i = 1:length(K_vector)
K_temp = K_vector(i);  
% Objective Function Pr.2
%fun = @(x) x(1).^2*K_temp/g  - (mu + 0.01*e)/(1-0.01*mu*e)  

fun = @(x) mu - x(1).^2*K_temp./g  
%C.1 (Bounds)
lb = 30; 
%ub = [30,35];
ub = 38;  % 60 < x2 < 80; mph
% There are no linear constraints, so set those arguments to |[]|. 
A = [];  b = []; % Linear In-equality Constraints
Aeq = []; beq = [];  % Linear Equality Constraints
%Initial Conditions
x0 = 30;  
%Constraints as an annoynomous function
nonlcon = @EqConstraint;
options = optimoptions('fmincon','Display','off');
Op(i,:) = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);  
end
fprintf('Pr. 2 Has finalized \n');
figure; plot(snew,Op(:,1),'linewidth',2)

xlabel('Segment Length s (m)'); ylabel('Optimized Speed (m/s)')
title('Segment Length vs Velocity Optimized'); grid on
%clc;
%cab(7)

% Nonlinear Constaints (Not bounds)
function [c,ceq] = EqConstraint(x)
global K_temp e g mu 
%Pr.2
% Nonlinear Inequality Constraints
c = x(1)^2*K_temp/g - (mu + 0.01*e)/(1-0.01*mu*e);
% Nonlinear Equality Constraints
ceq = [];
end



