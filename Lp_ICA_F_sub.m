function W = Lp_ICA_F_sub(X, D_sub, W_prev)
% Lp ICA - fast algorithm (O(n))
% X : PCA whitened data
% N : number of samples
% D_sub: number of sub-Gaussian sources
% D_sup: number of super-Gaussian sources

[D, N] = size(X);
if nargin < 2
    D_sub = D;
end

W = [];
if nargin == 3
    X = X - W_prev*(W_prev'*X);
    W = W_prev;
end

for ns = 1:D_sub
    % fast sub-gaussian algorithm
    % finding unmixing vector w by combinatorial optimization
    if ns > 1
        X_new = X - W*(W'*X);
    else
        X_new = X;
    end
    
    for i=1:N
        f(i) = norm(X_new(:,i));
    end
    
    [f_sort, ind] = sort(f,'descend');
    sgn = (round(rand(N,1))*2-1);
    sum_x = X_new(:,ind)*sgn;
    w = sum_x/norm(sum_x);
    a = X_new'*w;
    val = norm(a,4);
    
    for i=2:N
        sum_x1 = sum_x - 2*sgn(i)*X_new(:,ind(i));
        w = sum_x1/norm(sum_x1);
        a = X_new'*w;
        val1 = norm(a,4);
        if val1 <= val
            sum_x = sum_x1;
            val = val1;            
        end        
    end
    w = sum_x/norm(sum_x);    
    W = [W,w];    
end

end