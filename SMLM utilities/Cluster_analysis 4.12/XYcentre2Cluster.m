function [ X_cluster, Y_cluster] = XYcentre2Cluster( X_centre, Y_centre, N_loc_per_cluster, Cluster_size, Fov)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Display_info = 0;

N_clusters = length(X_centre);

X_cluster_all = zeros(N_clusters*N_loc_per_cluster, 1);
Y_cluster_all = zeros(N_clusters*N_loc_per_cluster, 1);

for i = 1:N_clusters 
    X_cluster_all((1+(i-1)*N_loc_per_cluster):(i*N_loc_per_cluster)) = X_centre(i) + (Cluster_size/2)*randn(N_loc_per_cluster,1);
    Y_cluster_all((1+(i-1)*N_loc_per_cluster):(i*N_loc_per_cluster)) = Y_centre(i) + (Cluster_size/2)*randn(N_loc_per_cluster,1);
end

X_cluster = X_cluster_all((X_cluster_all>0) & (Y_cluster_all>0) & (X_cluster_all<Fov) & (Y_cluster_all<Fov));
Y_cluster = Y_cluster_all((X_cluster_all>0) & (Y_cluster_all>0) & (X_cluster_all<Fov) & (Y_cluster_all<Fov));

if Display_info == 1
    disp(['Simulated density: ',num2str(size(X_cluster, 1)*1e6 / (Fov^2)),' particle / um^2']);
end

end

