function W = Lp_ICA_F(X, D_sub, D_sup, N)
% Lp ICA - fast algorithm (O(n))
% X : PCA whitened data
% N : number of samples
% D_sub: number of sub-Gaussian sources
% D_sup: number of super-Gaussian sources

D = D_sub+D_sup;
if D_sub > 0
    W1 = Lp_ICA_F_sub(X,D_sub);
    W = Lp_ICA_F_sup(X,D_sup,W1);
else
    W = Lp_ICA_F_sup(X,D_sup);
end


end