function W = Lp_ICA_F2(X, D_sub, D_sup, N)
% Lp ICA - fast algorithm (O(n))
% X : PCA whitened data
% N : number of samples
% D_sub: number of sub-Gaussian sources
% D_sup: number of super-Gaussian sources

D = D_sub+D_sup;

max_iter = 1000; Tol = 1E-20;
alpha = 0.1/N;

% probability model of the source
% :: f(s) = magnitude*exp(-scale*|s|^p)
scale =  [
    1.4123,
    0.5000,
    0.2281,
    %0.1142
    ]; % p=1:4

magnitude = [
    0.7062,
    0.3989,
    0.3421,
    %0.3207
    ]; % p=1:4

W = [];


for ns = 1:D_sub
    % fast sub-gaussian algorithm
    % finding unmixing vector w by combinatorial optimization
    if ns > 1
        X_new = X - W*(W'*X);
    else
        X_new = X;
    end
    
    for i=1:N
        x = X_new(:,i);
        if norm(x) < 1E-10
            f3(i) = 0;
        else
            w = x/norm(x);
            a = X'*w;
            f3(i) = norm(a,3); % 3-norm
        end
    end
    
    [f3_sort, ind3] = sort(f3,'descend');
    sgn = zeros(N,1);
    sgn(ind3(1),1) = 1;
    sum_x = X_new(:,ind3(1));
    
    for i=2:N-ns+1
        sum_x1 = sum_x + X_new(:,ind3(i));
        w_sub1 = sum_x1/norm(sum_x1);
        sum_x2 = sum_x - X_new(:,ind3(i));
        w_sub2 = sum_x2/norm(sum_x2);
        
        a1 = X'*w_sub1;
        a2 = X'*w_sub2;
        if norm(a1,3) < norm(a2,3)
            sum_x = sum_x1;
        else
            sum_x = sum_x2;
        end
    end
    w = sum_x/norm(sum_x);
    
    
    % concatenation of unmixing vector into the unmixing matrix
    W = [W,w];
    
end


XX = X;
for ns = 1:D_sup
    % fast super-gaussian algorithm
    
    % orthogonalization by Gram-Schmidt
    if ns > 1 || D_sub > 0
        XX = X - W*W'*X;
    end
    
    %norm of each column of X
    normX = sqrt(sum(XX.^2,1));
    %reject pre-selected index
    normX(normX<1e-10) = NaN;
    
    w = XX./ repmat(normX,D,1);
    a = X'*w;
    a = sum(abs(a),1);
    [~,ind] = min(a);
    w = w(:,ind);
    % concatenation of unmixing vector into the unmixing matrix
    W = [W,w];
        
end