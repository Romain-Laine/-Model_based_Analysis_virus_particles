function [ X1, Y1, X2_min, Y2_min] = AssociateCoordinates( X1, Y1, X2, Y2, R_search)
%ASSOCIATECOORDINATES Summary of this function goes here
%   Detailed explanation goes here

N_loc1 = size(X1,1);
N_loc2 = size(X2,1);

R_min = zeros(N_loc1,1);
X2_min = zeros(N_loc1,1);
Y2_min = zeros(N_loc1,1);
N_local = zeros(N_loc1,1);

for i = 1:N_loc1
    R = zeros(N_loc2,1);
    for j = 1:N_loc2
        R(j) = sqrt((X1(i)-X2(j))^2 + (Y1(i)-Y2(j))^2);
    end
    
    N_local(i) = sum(R < R_search); % Calculate the number of local minimum as defined by R_search
    R_min(i) = min(R);
    X2_min(i) = X2(R == min(R)); % Find the nearest localization
    Y2_min(i) = Y2(R == min(R));
end

X1(N_local ~= 1) = [];
Y1(N_local ~= 1) = [];

X2_min(N_local ~= 1) = [];
Y2_min(N_local ~= 1) = [];

% Display
figure('Color','white','name','Histogram of R_min','Units','normalized','OuterPosition',[0.2 0.2 0.6 0.5]);
subplot(1,2,1)
plot(X1,Y1,'+')
hold on
plot(X2_min,Y2_min,'r+')
axis equal
legend('Red channel','Green channel');

subplot(1,2,2)
hist(R_min,0:10:500)
xlim([0 500])
xlabel 'R_{offset} (nm)'
title 'Histogram of chromatic offset'


end

