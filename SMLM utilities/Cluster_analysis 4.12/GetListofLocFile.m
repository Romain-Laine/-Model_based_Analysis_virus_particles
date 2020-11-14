function [ PathName, FileNames ] = GetListofLocFile( DefaultPath, File_token, area_token )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

PathName = uigetdir(DefaultPath , 'Choose directory containing the localization files (.txt)...');
FileList = dir([PathName,'\',area_token,'*',File_token,'.txt']);

N_files = length(FileList);
FileNames = cell(N_files,1);

for i = 1:N_files
    FileNames{i} = [PathName, '\', FileList(i).name];
end


end

