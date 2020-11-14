function [ DirPath, FileNames_2Cdata ] = get2ColorLocFiles(DefaultPathName, RC_token, GC_token )
%GET2COLORLOCFILES Summary of this function goes here
%   Detailed explanation goes here

DirPath = uigetdir(DefaultPathName, 'Choose a folder containing the localization files...');
FileList = dir([DirPath,'\*.txt']);

FileNames = cell(size(FileList,1),1);
for i = 1:size(FileList,1)
    FileNames{i} = FileList(i).name;
end

Found = 1;
area_number = 0;
while Found;
    area_number = area_number + 1;
    Found_RC = (ismember(['a',num2str(area_number),RC_token,'.txt'],FileNames));
    Found_GC = (ismember(['a',num2str(area_number),GC_token,'.txt'],FileNames));
    Found = Found_RC && Found_GC;
end

area_number = area_number - 1;

FileNames_2Cdata = cell(area_number,1);
for i = 1:area_number
    FileNames_2Cdata{i,1} = ['a',num2str(i),RC_token,'.txt'];
    FileNames_2Cdata{i,2} = ['a',num2str(i),GC_token,'.txt'];
end


end

