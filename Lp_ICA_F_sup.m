function W = Lp_ICA_F_sup(X, D_sup, W_prev)
% Lp ICA - fast algorithm (O(n))
% X : PCA whitened data
% N : number of samples
% D_sub: number of sub-Gaussian sources
% D_sup: number of super-Gaussian sources

[D, N] = size(X);
if nargin < 2
    D_sup = D;
end

W = [];
if nargin == 3
    X = X - W_prev*(W_prev'*X);
    W = W_prev;
end

for ns = 1:D_sup
    % fast sub-gaussian algorithm
    % finding unmixing vector w by combinatorial optimization
    if ns > 1
        X_new = X - W*(W'*X);
    else
        X_new = X;
    end
   
    for i=1:N
        x = X_new(:,i);
        if norm(x) < 1E-10  % if x == 0, we do not care the sample
            f(i) = 1E10*N; 
            %f3(i) = norm(x);
        else
            w = x/norm(x);
            a = X_new'*w;
            f(i) = norm(a,1);  % 1-norm
        end
    end
    
    [f_min, ind] = min(f);
    x = X_new(:,ind);
    w = x/norm(x);    % optimal projection vector for super-G source
    
    W = [W,w];    
end

end