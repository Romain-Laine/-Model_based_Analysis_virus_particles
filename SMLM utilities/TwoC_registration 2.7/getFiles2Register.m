function [ DirPath, FileNames ] = getFiles2Register(DefaultPathName, GC_token )
%GET2COLORLOCFILES Summary of this function goes here
%   Detailed explanation goes here

DirPath = uigetdir(DefaultPathName, 'Choose a folder containing the localization files...');
FileList = dir([DirPath,'\*',GC_token ,'.txt']);

FileNames = cell(size(FileList,1),1);
for i = 1:size(FileList,1)
    FileNames{i} = FileList(i).name;
end

end

