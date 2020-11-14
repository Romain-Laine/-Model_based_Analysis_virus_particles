function [ ParticleInfo, FolderName ] = Multiple_COM( PathName )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

FolderName = uigetdir(PathName,'Select a folder: ');
List = dir([FolderName,'\ParticleInfo_*.txt']);

nObject = 0;
ParticleInfo = [0 0 0 0 0];
for i = 1:max(size(List))
    ReadInfo = dlmread([FolderName,'\',List(i).name],'\t',1,0);
    nObject = cat(1,nObject,ReadInfo(:,5)+max(nObject));
    ParticleInfo = cat(1,ParticleInfo,ReadInfo);
end

ParticleInfo(:,5) = nObject;
ParticleInfo(1,:) = [];

end

