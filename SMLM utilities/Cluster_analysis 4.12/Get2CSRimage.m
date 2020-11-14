function [ PathName, FileNames_RC,FileNames_GC ] = Get2CSRimage( DefaultPath, RC_token, GC_token, area_token )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

PathName = uigetdir(DefaultPath , 'Choose directory containing the png STORM images...');
FileList_RC = dir([PathName,'\',area_token,'*',RC_token,'.png']);
FileList_GC = dir([PathName,'\',area_token,'*',GC_token,'.png']);

N_files_RC = length(FileList_RC);
N_files_GC = length(FileList_GC);

if N_files_RC == N_files_GC
    FileNames_RC = cell(N_files_GC,1);
    FileNames_GC = cell(N_files_GC,1);
    
    for i = 1:N_files_RC
        FileNames_RC{i} = [PathName, '\', FileList_RC(i).name];
        FileNames_GC{i} = [PathName, '\', FileList_GC(i).name];
    end
end



end

