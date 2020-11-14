% Romain Laine - Laser Analytics Group
% rfl30@cam.ac.uk
% This script will use the tform trafo saved from TwoC_reg.m to register localization data from the GC.
% The localization data are then saved in a loc file with the correct format to be used in rapidSTORM replay.
% The format is x y t I
% -------------------------------------------------------------------------

clear all
close all
clc

% Header from rapidSTORM
Header = '# <localizations insequence="true" repetitions="variable"><field identifier="Position-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position in sample space in X" unit="nanometer" min="0 m" max="4.0086e-005 m" /><field identifier="Position-1-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="position in sample space in Y" unit="nanometer" min="0 m" max="4.0086e-005 m" /><field identifier="ImageNumber-0-0" syntax="integer" semantic="frame number" unit="frame" min="0 fr" /><field identifier="Amplitude-0-0" syntax="floating point with . for decimals and optional scientific e-notation" semantic="emission strength" unit="A/D count" min="10000 ADC" /></localizations>';
DefaultPathName = 'C:\Users\rfl30\DATA raw\dSTORM v1 data\';

GC_token = '_BC';

% Load the trafo
[TrafoFileName,TrafoPathName] = uigetfile([DefaultPathName, '\*.mat'],'Choose the transformation file');
load([TrafoPathName,TrafoFileName]); % load and save as tform in the workspace

% Find the data to register in a user defined folder
[ DirPath, FileNames ] = getFiles2Register(TrafoPathName, GC_token );
disp('Files found in the folder:');
disp(FileNames);
Save_button = questdlg('Save registered localization file?','Save localization file','Yes','No','No');


for i = 1:size(FileNames,1)
    disp('------------------------------');
    
    % Load the data to correct for abberation (Green Channel)
    disp([DirPath,'\',FileNames{i}]);
    LocInfo = Read_LocFile([DirPath,'\',FileNames{i}]);
    X_d = LocInfo(:,1);
    Y_d = LocInfo(:,2);
    
    % Apply the obtained trafo
    [U,V] = tforminv(tform,X_d,Y_d);
    
    % Display
    figure('Color','white','name',['Scatter plot ( ',FileNames{i},' )']);
    plot(X_d,Y_d,'+')
    hold on
    plot(U,V,'r+')
    axis equal
    legend('Raw','Registered');
    
    
    %% Save the registered localization file
    if strcmp(Save_button,'Yes')
        N_line = size(U,1);
        % new file name
        NewFileName = [regexprep(FileNames{i},'.txt',''),'_registered.txt'];
        
        % Write on the file with the correct format for rapidSTORM
        fid = fopen([DirPath,'\',NewFileName],'wt');
        fprintf( fid, Header);
        for j = 1:N_line
            fprintf( fid, '%.2f %.2f %i %.2f\n', U(j), V(j), LocInfo(j,3), LocInfo(j,4));
        end
        
        disp('Saved as...');
        disp([DirPath,'\',NewFileName]);
        
        fclose(fid);
    end
    
    
end

disp('------------------------------');
disp('All done.')





