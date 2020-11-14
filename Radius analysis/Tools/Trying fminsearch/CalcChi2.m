function [ Chi2 ] = CalcChi2( rt,r_hist,n )

r = rt(1);
t = rt(2);

N = 1000000;
xy = MC_Sim_3DShell(N, r, t);
r_sim = sqrt(xy(:,1).^2 + xy(:,2).^2);
n_sim = hist(r_sim,r_hist);
n_sim = sum(n)*n_sim/(sum(n_sim));

W = n;
W(W == 0) = 1;
Chi2 = (1/(size(r_hist,2) - 2 - 1))*sum((n - n_sim).^2./W);

end

