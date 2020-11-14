function [ ParticleInfo, FolderName ] = Multiple_COM( PathName )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

FolderName = uigetdir(PathName,'Select a folder: ');
List = dir([FolderName,'\ParticleInfo_*.txt']);

nObject = 0;  % Object number
ParticleInfo = [0 0 0 0 0];

for i = 1:max(size(List))
    ReadInfo = dlmread([FolderName,'\',List(i).name],'\t',1,0);
    nObject = cat(1,nObject,ReadInfo(:,5)+max(nObject));
    ParticleInfo = cat(1,ParticleInfo,ReadInfo);
end

% Append the column with particle numbers
ParticleInfo(:,5) = nObject;

% Get rid of the first line that's only got zeros in
ParticleInfo(1,:) = [];

ADCs = ParticleInfo(:,4);

disp('-----------------------');
disp('Multiple particle ADC statistics:');
disp(['Mean: ',num2str(round(mean(ADCs))),' ADC']);

% Round up to the 10
roundADCs = 10*round(ADCs/10);
disp(['Mode: ',num2str(mode(roundADCs)),' ADC']);
disp(['Median: ',num2str(median(roundADCs)),' ADC']);

% ADC_hist = linspace(min(ADCs),10*median(roundADCs),round(sqrt(size(ADCs,1))));
% figure('color','white','name','ADC histogram from multiple particles');
% hist(ADCs,ADC_hist);
% xlabel 'ADC'
% ylabel 'Occurences'
% title(['Average ADC: ',num2str(round(mean(ADCs))),' ADC'])

end

