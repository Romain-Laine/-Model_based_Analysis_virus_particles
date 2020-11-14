clear all
close all
clc

ProtStructure = 'Shell';
% ProtStructure = 'Ring';

LinkerStructure = 'Shell';
% LinkerStructure = 'Ring';

Linker_size = 10;       % nm
Linker_thickness = 2;   % nm
LocError = 20;          % nm

N_sim = 1000000;

N_param = 5;
% r = linspace(60,100,N_param);
% t = 20*ones(1,N_param);

r = 80*ones(1,N_param);
t = linspace(10,30,N_param);


r_hist = 0:1:160;
Sim_param = cell(5,1);
Sim_param{1} = ProtStructure;
Sim_param{2} = LinkerStructure;
Sim_param{3} = Linker_size;        % nm
Sim_param{4} = Linker_thickness;   % nm
Sim_param{5} = LocError;           % nm

n_opt = zeros(size(r_hist,2),N_param);

h = waitbar(0,'Please wait...');
for i = 1:N_param
    n_opt(:,i) = VirusSim(Sim_param, r(i), t(i), N_sim, r_hist);
    waitbar(i/N_param,h);
end

close(h);

%%
figure('Color','white');
plot(r_hist,n_opt)
xlabel 'r (nm)'
ylabel 'Occurences'
grid on

disp(cat(1,r,t))

r_model = 80;
t_model = 10;
Shell_model = Shell_Analytical_Model(r_hist,r_model,t_model,LocError,Linker_size, Linker_thickness);
Shell_model = Shell_model*sum(n_opt(:,1));

figure('Color','white');
plot(r_hist,n_opt)
xlabel 'r (nm)'
ylabel 'Occurences'
grid on
hold on
plot(r_hist,Shell_model,'g--')



