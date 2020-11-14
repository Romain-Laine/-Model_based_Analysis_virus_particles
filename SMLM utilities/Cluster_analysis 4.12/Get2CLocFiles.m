function [ FileNames_RC,FileNames_GC ] = Get2CLocFiles( DefaultPath )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

RC_token = '_647_stack';
GC_token = '_568_stack_registered';

PathName = uigetdir(DefaultPath , 'Choose directory conataing the localization files...');
FileList_RC = dir([PathName,'\area*',RC_token,'.txt']);
FileList_GC = dir([PathName,'\area*',GC_token,'.txt']);

N_files_RC = size(FileList_RC,1);
N_files_GC = size(FileList_GC,1);


if N_files_RC == N_files_GC
    FileNames_RC = cell(N_files_GC,1);
    FileNames_GC = cell(N_files_GC,1);
    
    for i = 1:N_files_RC
        FileNames_RC{i} = [PathName, '\', FileList_RC(i).name];
        FileNames_GC{i} = [PathName, '\', FileList_GC(i).name];
    end
end



end

