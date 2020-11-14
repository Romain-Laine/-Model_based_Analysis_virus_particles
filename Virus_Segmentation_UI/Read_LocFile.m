function LocInfo = Read_LocFile( FileName )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Header
% semantic="position in sample space in x dimension" unit="nanometer" min="0 m" max="4.0392e-005 m" />
% semantic="position in sample space in y dimension" unit="nanometer" min="0 m" max="4.0392e-005 m" />
% semantic="frame number" unit="frame" min="0 fr" /
% semantic="emission strength" unit="A/D count" />
% semantic="two kernel improvement" unit="dimensionless" />
% semantic="fit residue chi square value" unit="dimensionless"

block_size = 25000;
format = '%f %f %f %f %*[^\n]';   % %*[^\n] ignores the rest of the line

tic
file_id = fopen(FileName);
segarray = textscan(file_id, format,block_size, 'delimiter',' ','CommentStyle','#');
LocInfo = cell2mat(segarray);


while ~feof(file_id)
    segarray = textscan(file_id, format, block_size, 'delimiter',' ','CommentStyle','#');
    LocInfo = cat(1,LocInfo,cell2mat(segarray));
end

fclose(file_id);
toc

disp(['Number of localizations: ',num2str(size(LocInfo,1))]);



end

