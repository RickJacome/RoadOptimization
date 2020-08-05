%   Based off the code from 
%   Are Mjaavatten (2019). Curvature of a 2D or 3D curve
%   (https://www.mathworks.com/matlabcentral/
%   fileexchange/69452-curvature-of-a-2d-or-3d-curve),
%   MATLAB Central File Exchange. 
function [L,R,kappa] = curvature(X)
% Radius of curvature and curvature vector for 2D or 3D curve
%  [L,R,Kappa] = curvature(X)
%   X:   2 or 3 column array of x, y (and possibly z) coordiates
%   L:   Cumulative arc length
%   R:   Radius of curvature
%   k:   Curvature vector

  n = size(X,1);
  dims = size(X,2); %Check Dimension Size
  if dims == 2
    X = [X,zeros(n,1)];  % Do all calculations in 3D
  end
  L = zeros(n,1);
  R = NaN(n,1);  %NaN Vector Size of N
  kappa = NaN(n,3); 
  for i = 2:n-1
      %The input for the Circumcenter function 
      %is X current (BOTH coordinates of 1 point), X previous (BOTH
      % coordinates of second point, and X future (both coordinates of
      % third point)
      % Outputs are Radius R and Curvature k
    [R(i),~,kappa(i,:)] = circumcenter(X(i,:)',X(i-1,:)',X(i+1,:)');
    L(i) = L(i-1)+norm(X(i,:)-X(i-1,:));
  end
  i = n;
  L(i) = L(i-1)+norm(X(i,:)-X(i-1,:));
  if dims == 2
    kappa = kappa(:,1:2);
  end
end

