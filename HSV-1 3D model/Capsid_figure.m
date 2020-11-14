clear all
close all
clc

t   = (1+sqrt(5))/2;
tau = t/sqrt(1+t^2);
one = 1/sqrt(1+t^2);

p = [
    +tau  +one  +0     % ZA
    -tau  +one  +0     % ZB
    -tau  -one  +0     % ZC
    +tau  -one  +0     % ZD
    +one  +0    +tau   % YA
    +one  +0    -tau   % YB
    -one  +0    -tau   % YC
    -one  +0    +tau   % YD
    +0    +tau  +one   % XA
    +0    -tau  +one   % XB
    +0    -tau  -one   % XC
    +0    +tau  -one]; % XD

Theta_x = -10*pi/180;
Rx = [1 0 0; 0 cos(Theta_x) -sin(Theta_x); 0 sin(Theta_x) cos(Theta_x)];

Theta_y = -10*pi/180;
Ry = [cos(Theta_y) 0 sin(Theta_y); 0 1 0; -sin(Theta_y) 0 cos(Theta_y)];


R = Rx*Ry;

for i = 1:12
    p_r = R*p(i,:)';
    p(i,:) = p_r';
end


t = [
    5  8  9
    5 10  8
    6 12  7
    6  7 11
    1  4  5
    1  6  4
    3  2  8
    3  7  2
    9 12  1
    9  2 12
    10  4 11
    10 11  3
    9  1  5
    12  6  1
    5  4 10
    6 11  4
    8  2  9
    7 12  2
    8 10  3
    7  3 11 ];


figure(1);
axis equal
hold on
trisurf(t,p(:,1),p(:,2),p(:,3));
axis vis3d
view(3);


Max_xy = 1.2;
Step = 0.005;

gv = -Max_xy:Step:Max_xy;  % nm
[X,Y,Z] = meshgrid(gv);

Cube = ones(size(X));
n_vectors = zeros(20,3);

for s = 1:20
    n = mean([p(t(s,1),:); p(t(s,2),:); p(t(s,3),:)]);
    n_vectors(s,:)  = n/norm(n);
end

figure;
for s = 1:20
    quiver3(0,0,0,n_vectors(s,1),n_vectors(s,2),n_vectors(s,3));
    hold on
end
axis equal

%%

tic
disp('Cutting off planes...')
h_wait = waitbar(0,'Please wait...');

for s = 1:20
    waitbar(s/20);
    Xn = X - n_vectors(s,1);
    Yn = Y - n_vectors(s,2);
    Zn = Z - n_vectors(s,3);
    
    D = Xn*n_vectors(s,1) + Yn*n_vectors(s,2) + Zn*n_vectors(s,3);
    Cube(D > 0) = 0;
    
end

close(h_wait);
toc


%%

disp('Smoothing...');
tic
Cube = smooth3(Cube,'gaussian',11,4);
% Cube(Y < 0 & Z > 0 & X > 0) = 0;   % this line cuts off a bit of the
% capsid
toc

% Display everything

figure('Color','white');

disp('Isosurface...');
tic
isoval = 0.5;
hisoxy = patch(isosurface(Cube,isoval),...
    'FaceColor',[1 0.5 0],...
    'EdgeColor','none');
isonormals(Cube,hisoxy)
toc

set(gcf,'Renderer','zbuffer'); lighting phong
set(hisoxy,'SpecularColorReflectance',0.2,'SpecularExponent',50)

lightangle(30,30);
view(24, 12)
axis off
grid off
axis equal