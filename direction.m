function [O1,O2] = direction(K2)
    %Orthogonal Shift Function 
     O1 = atand(K2(:,2)./K2(:,1));  %Angle in Degrees
     O2 = O1 + 90;   % 90 Degrees
end