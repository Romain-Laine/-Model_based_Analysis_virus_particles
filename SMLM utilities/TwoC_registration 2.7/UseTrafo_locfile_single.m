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
DefaultPathName = 'E:\STORM data\2014_10_27_2Cregistration\Averaged data';

% Load the trafo
[FileName,PathName] = uigetfile([DefaultPathName, '\*.mat'],'Choose the transformation file');
load([PathName,FileName]);

% Load the data to correct for abberation (Green Channel)
[FileName, PathName] = uigetfile([DefaultPathName,'\*.txt'],'Select green channel coordinates dataset');
disp([PathName,FileName]);

LocInfo = Read_LocFile([PathName,FileName]);
X_d = LocInfo(:,1);
Y_d = LocInfo(:,2);

% Apply the obtained trafo
[U,V] = tforminv(tform,X_d,Y_d);

% Display
figure('Color','white','name','Scatter plot (raw)');
plot(X_d,Y_d,'+')
hold on
plot(U,V,'r+')
axis equal


%% Save the registered localization file
button = questdlg('Save registered localization file?','Save localization file','Yes','No','No');
if strcmp(button,'Yes')
    N_line = size(U,1);
    % new file name
    NewFileName = [regexprep(FileName,'.txt',''),'_registered.txt'];
    
    % Write on the file with the correct format for rapidSTORM
    fid = fopen([PathName,NewFileName],'wt');
    fprintf( fid, Header);
    for i = 1:N_line
        fprintf( fid, '%.2f %.2f %i %.2f\n', U(i), V(i), LocInfo(i,3), LocInfo(i,4));
    end
    
    fclose(fid);
end






