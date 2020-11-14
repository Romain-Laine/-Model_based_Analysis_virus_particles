clear all
close all
clc

R_search = 200;
DefaultPathName = 'E:\STORM data\2014_10_27_2Cregistration\Averaged data';

% Load the trafo
[FileName,PathName] = uigetfile([DefaultPathName, '\*.mat'],'Choose the transformation file');
load([PathName,FileName]);

% Load the data to correct for abberation (Green Channel)
[FileName, PathName] = uigetfile([DefaultPathName,'\*.txt'],'Select red channel coordinates dataset');
disp([PathName,FileName]);

ReadInfo_RC = dlmread([PathName,FileName],' ',1,0);
X = ReadInfo_RC(:,1);
Y = ReadInfo_RC(:,2);
disp(['Number of localizations in GC: ',num2str(size(X,1))]);

[FileName, PathName] = uigetfile([DefaultPathName,'\*.txt'],'Select green channel coordinates dataset');
disp([PathName,FileName]);

ReadInfo_GC = dlmread([PathName,FileName],' ',1,0);
X_d = ReadInfo_GC(:,1);
Y_d = ReadInfo_GC(:,2);
disp(['Number of localizations in RC: ',num2str(size(X_d,1))]);

disp('-----------------');
[X, Y, X_d, Y_d] = AssociateCoordinates(X, Y, X_d, Y_d, R_search);
disp(['Number of localizations associated: ',num2str(size(X,1))]);

% Display
figure('Color','white','name','Scatter plot (raw)');
plot(X,Y,'+')
hold on
plot(X_d,Y_d,'r+')
axis equal


% Apply the obtained trafo 
[U,V] = tforminv(tform,X_d,Y_d);

figure('Color','white','name','Post-registration scatter plot');
plot(X,Y,'+')
hold on
plot(U,V,'r+')
axis equal

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










