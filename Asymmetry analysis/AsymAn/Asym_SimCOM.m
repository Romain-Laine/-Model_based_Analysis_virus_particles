clear all
close all
clc

N_particles = 2;
LocError = 10;
r_sim = 15;
N_repeat = 50000;

D = zeros(1,N_repeat);
x0_com = zeros(1,N_repeat);
y0_com = zeros(1,N_repeat);

tic
for i = 1:N_repeat
    % Generate a dataset with Thickness = 0
    xy_sim = MC_Sim_3DShell(N_particles,r_sim,0) + MC_Sim_2DLocError(N_particles,LocError);
    x0_com(i) = mean(xy_sim(:,1));
    y0_com(i) = mean(xy_sim(:,2));
end

% figure('Color','white');
% hist(x0_com,100);


%%

xy0 = cat(2,x0_com,y0_com);
figure('Color','white');
hist(xy0,100);

disp('-------------------------');
disp(['Number of particles: ',num2str(N_particles)]);
disp(['Localization error : ',num2str(LocError),' nm']);
disp('-------------------------');
disp(['Mean = ',num2str(mean(xy0)), ' nm']);
disp(['STD  = ',num2str(std(xy0)), ' nm']);


toc