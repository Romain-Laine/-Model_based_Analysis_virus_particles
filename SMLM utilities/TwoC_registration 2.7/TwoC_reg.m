% RFL 2014-10-28 rfl30@cam.ac.uk
% This code can simulate some bead data or take some measured ones to generate a trafo.
% If simulated: random or regular pattern of beads can be generated.
% If measured: the localization coordinates from the two channels are associated.
% This association is performed by looking at how many neighbours are present within the circle defined by R_search

clear all
close all
clc

% Bead_type = 'Regular';
% Bead_type = 'Random';
Bead_type = 'Measured';

% ------------ For simulated data ----------------
n_points_reg = 50;

LocError = 10;
n_points = 20;
Field_size = 40000; % nm (for simulated data)

SimTrafo_type = 'Polynomial r -> r + a*r + b*r^2';
Defo_centre = [20000, 20000]; % nm
Defo_Amp_lin = 0.005;
Defo_Amp_quad = 0.000000;
xy_shift = [0, 0]; % nm

Trafo_param = cell(5,1);
Trafo_param{1} = SimTrafo_type;
Trafo_param{2} = Defo_centre;
Trafo_param{3} = Defo_Amp_lin;
Trafo_param{4} = Defo_Amp_quad;
Trafo_param{5} = xy_shift;

% ------------ For measured data ----------------
DefaultPathName = 'C:\Users\rfl30\DATA raw\dSTORM v1 data\2015_10_01 Micelles\Beads registration\';
R_search = 200; % in nm

RC_token = '_RC_AVG';     % file found will be: a1_RC_AVG.txt, a2_RC_AVG.txt etc. 
GC_token = '_BC_AVG';

% ------------ For analysing the trafo ----------------

% Trafo_type = 'lwm';
Trafo_type = 'polynomial';

N_lwm = 12;
polyn_order = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtaining the positions of the beads in the image

if strcmp(Bead_type,'Regular')
    xgv = linspace(0,Field_size,n_points_reg);
    [X,Y] = meshgrid(xgv,xgv);
    X = reshape(X,[size(X,1)*size(X,2),1]);
    Y = reshape(Y,[size(Y,1)*size(Y,2),1]);
    [ X_d, Y_d ] = SimTrafo( X,Y, Trafo_param);
    
    % apply some localization error
    [ xy ] = MC_Sim_2DLocError(size(X,1),LocError);
    X = X + xy(:,1);
    Y = Y + xy(:,2);
    [ xy ] = MC_Sim_2DLocError(size(X_d,1),LocError);
    X_d = X_d + xy(:,1);
    Y_d = Y_d + xy(:,2);
    
    
elseif strcmp(Bead_type,'Random')
    X = rand(n_points_reg,1)*Field_size;
    Y = rand(n_points_reg,1)*Field_size;
    [ X_d, Y_d ] = SimTrafo( X,Y, Trafo_param);
    
    % apply some localization error
    [ xy ] = MC_Sim_2DLocError(size(X,1),LocError);
    X = X + xy(:,1);
    Y = Y + xy(:,2);
    [ xy ] = MC_Sim_2DLocError(size(X_d,1),LocError);
    X_d = X_d + xy(:,1);
    Y_d = Y_d + xy(:,2);
    
    
    % From measured data
elseif strcmp(Bead_type,'Measured')
    
    [ PathName, FileNames] = get2ColorLocFiles(DefaultPathName, RC_token, GC_token );
    if isempty(FileNames)
        disp('No files complementary files found...');
        return
    end
    
    disp('Files found in folder:');
    disp(FileNames);
    
    N_file = size(FileNames,1);
    
    X = [];
    Y = [];
    X_d = [];
    Y_d = [];
    
    for i = 1:N_file
        disp('-----------------');
        disp([PathName,'\',FileNames{i,1}]);
        ReadInfo_RC = dlmread([PathName,'\',FileNames{i,1}],' ',1,0);
        Xi = ReadInfo_RC(:,1);
        Yi = ReadInfo_RC(:,2);
        disp(['Number of localizations in GC: ',num2str(size(Xi,1))]);
        
        disp([PathName,'\',FileNames{i,2}]);
        ReadInfo_GC = dlmread([PathName,'\',FileNames{i,2}],' ',1,0);
        X_di = ReadInfo_GC(:,1);
        Y_di = ReadInfo_GC(:,2);
        disp(['Number of localizations in RC: ',num2str(size(X_di,1))]);
        
        [Xi, Yi, X_di, Y_di] = AssociateCoordinates(Xi, Yi, X_di, Y_di, R_search);
        set(gcf,'name',['Scatter plot (associated) (',FileNames{i,1},' and ',FileNames{i,2},')']);
        disp(['Number of localizations associated: ',num2str(size(Xi,1))]);
        
        X = cat(1,X,Xi);
        Y = cat(1,Y,Yi);
        X_d = cat(1,X_d,X_di);
        Y_d = cat(1,Y_d,Y_di);
        
    end
    
end


% Estimate transformation -------------------------------------------------

disp(' ');
disp('---------------------------- Results ---------------------------------');
disp(['Trafo type: ',Trafo_type]);
if strcmp(Trafo_type,'polynomial')
    param = polyn_order;
elseif strcmp(Trafo_type,'lwm')
    param = N_lwm;
end

tform = cp2tform([X Y],[X_d Y_d],Trafo_type,param);

% Apply the obtained trafo
[U,V] = tforminv(tform,X_d,Y_d);

% display
figure('Color','white','name','Scatter plot pre- and post-registration','Units','normalized','Outerposition',[0.1 0.1 0.8 0.6]);
subplot(1,2,1)
plot(X,Y,'+')
hold on
plot(X_d,Y_d,'r+')
axis equal
title 'Pre-registration'
legend('Red channel','Green channel');

subplot(1,2,2)
plot(X,Y,'+')
hold on
plot(U,V,'r+')
axis equal
title 'Post-registration'
legend('Red channel','Green channel');

% Calculate the TRE
PreRegOffset = sqrt((X-X_d).^2 + (Y-Y_d).^2);
TRE_pre = mean(PreRegOffset);
disp(['Pre-reg TRE = ',num2str(TRE_pre,'%6.1f'),' nm']);

PostRegOffset = sqrt((X-U).^2 + (Y-V).^2);
TRE = mean(PostRegOffset);
disp(['Post-reg TRE = ',num2str(TRE,'%6.1f'),' nm']);

figure('Color','white','name','Histogram of Post-registration offset');
hist(PostRegOffset,0:2:50)
xlim([0 50])
xlabel 'R_{offset} (nm)'


%% Testing trafo -----------------------------------------------------------

xgv = linspace(0,Field_size,n_points);
xgv(1) = [];
xgv(end) = [];

[X,Y] = meshgrid(xgv,xgv);
X = reshape(X,[size(X,1)*size(X,2),1]);
Y = reshape(Y,[size(Y,1)*size(Y,2),1]);

if strcmp(Bead_type,'Measured')
    
    % Deformation -------------------------------------------------------------
    [X_dinv,Y_dinv] = tforminv(tform,X,Y);
    figure('Color','white','name','Displaying trafo','Units','normalized','Outerposition',[0.1 0.1 0.8 0.6]);
    subplot(1,2,1)
    plot(X,Y,'+')
    hold on
    plot(X_dinv,Y_dinv,'r+')
    axis equal
    
    subplot(1,2,2)
    quiver(X,Y,X -X_dinv, Y-Y_dinv)
    axis equal
    
    
else
    disp('Testing simulated trafo...')
    
    % Deformation -------------------------------------------------------------
    [ X_d, Y_d ] = SimTrafo( X,Y, Trafo_param);
    
    [U,V] = tforminv(tform,X_d,Y_d);
    figure('Color','white','name','Testing trafo - Registered coordinates');
    plot(X,Y,'+')
    hold on
    plot(U,V,'r+')
    axis equal
    xlim([0 Field_size])
    ylim([0 Field_size])
    
    figure('Color','white','name','Displacement map');
    quiver(X,Y,X -X_d, Y-Y_d)
    axis equal
    
    
    % Display offset
    
    X = reshape(X,[n_points,n_points]);
    Y = reshape(Y,[n_points,n_points]);
    U = reshape(U,[n_points,n_points]);
    V = reshape(V,[n_points,n_points]);
    R = sqrt((U-X).^2 + (V-Y).^2); % in nm
    
    %%
    figure('Color','white','name','Testing trafo - TRE image');
    imshow(R,[])
    title(['# point in registration: ',num2str(n_points_reg^2),' - Max offset: ',num2str(max(R(:)),'%6.1f'),' nm - FRE = ',num2str(mean(R(:)),'%6.1f'),' nm'])
    
    figure('Color','white','name','Testing trafo - Line profile through centre of TRE image');
    plot(xgv,R(round(size(R,1)/2),:))
    xlabel 'um'
    ylabel 'Offset (nm)'
    title(['# point in registration: ',num2str(n_points_reg^2),' - Max offset: ',num2str(max(R(:)),'%6.1f'),' nm - FRE = ',num2str(mean(R(:)),'%6.1f'),' nm'])
end

%%
button = questdlg('Save transformation?','Save transformation','Yes','No','No');
if strcmp(button,'Yes')
    save([PathName,'\Trafo_tform.mat'],'tform');
    disp(['Transformation saved as: ',PathName,'\Trafo_tform.mat']);
end






