function [ Hist_sim_out ] = VirusSim(Sim_param, radVP, ThicknessVP, N_loc_total, r_hist)
% VirusSim is a function based on AlignedVirusSIM m-file to simulate a
% shell or ring of virus protein and their localizations


ProtStructure = Sim_param{1};
% ProtStructure = 'Ring' or 'Shell';
LinkerStructure = Sim_param{2};
% LinkerStructure = 'Ring' or 'Shell';

Linker_size = Sim_param{3};        % nm
Linker_thickness = Sim_param{4};   % nm
LocError = Sim_param{5};           % nm

Chunk_size = 100000000;
n_chunks = ceil(N_loc_total/Chunk_size);

N_loc = min(N_loc_total, Chunk_size);
%------------------------------------------------------------------------------------------------------------------------%
% Generate the protein localization
Hist_sim = zeros(n_chunks, size(r_hist,2));

for i = 1:n_chunks
    if strcmp(ProtStructure,'Shell')
        Prot_loc = MC_Sim_3DShell(N_loc, radVP, ThicknessVP);
    elseif strcmp(ProtStructure,'Ring')
        Prot_loc = MC_Sim_2DRing(N_loc, radVP, ThicknessVP);
    end
    
    % Generate the fluorophore localization
    if strcmp(LinkerStructure,'Shell')
        Fluo_loc = Prot_loc + MC_Sim_3DShell(N_loc, Linker_size, Linker_thickness);
    elseif strcmp(LinkerStructure,'Ring')
        Fluo_loc = Prot_loc + MC_Sim_2DRing(N_loc, Linker_size, Linker_thickness);
    end
    
    % Generate the error on the localization due to precision
    xy_loc_error = MC_Sim_2DLocError(N_loc,LocError);
    
    % Sum it all up to get the localization
    xy = Fluo_loc + xy_loc_error;
    
    Hist_sim(i,:) = hist(sqrt((xy(:,1)).^2 + (xy(:,2)).^2),r_hist);
end

if n_chunks == 1
    Hist_sim_out = Hist_sim;
elseif n_chunks > 1
    Hist_sim_out = mean(Hist_sim);
end
    
    
end

