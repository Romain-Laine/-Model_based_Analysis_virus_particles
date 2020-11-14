function [ ParticleInfo, FolderName ] = Multiple_load_loc_file( PathName, Token_Name )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

FolderName = uigetdir(PathName,'Select a folder: ');
List = dir([FolderName,'\a*', Token_Name, '.txt']);
disp(['Number of files: ', num2str(size(List,1))]);

nObject = 0;
ParticleInfo = [0 0 0 0 0];
for i = 1:size(List,1)
    i
    ReadInfo = dlmread([FolderName,'\',List(i).name],'\t',1,0);
    size(ReadInfo)
    nObject = cat(1,nObject,ReadInfo(:,5)+max(nObject));
    ParticleInfo = cat(1,ParticleInfo,ReadInfo);
end

ParticleInfo(:,5) = nObject;
ParticleInfo(1,:) = [];

end

