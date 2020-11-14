function [ X,Y ] = Generate_clustered_XY( N_clusters, N_loc_clusters, Cluster_size, Fov )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x_clusters = Fov*rand(N_clusters,1);
y_clusters = Fov*rand(N_clusters,1);

x_all = zeros(0,1);
y_all = zeros(0,1);
for i = 1:N_clusters
    xi = x_clusters(i) + Cluster_size/2*randn(N_loc_clusters,1);
    yi = y_clusters(i) + Cluster_size/2*randn(N_loc_clusters,1);
    x_all = cat(1,x_all,xi);
    y_all = cat(1,y_all,yi);
end

X = x_all((x_all>0) & (y_all>0) & (x_all<Fov) & (y_all<Fov));
Y = y_all((x_all>0) & (y_all>0) & (x_all<Fov) & (y_all<Fov));
disp(['Simulated density: ',num2str(size(X,1)*1e6/(Fov^2)),' particle / um^2']);

end

