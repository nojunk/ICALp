function W = ICA_PCA_L1(X,D,N)
% ICA_PCA_L1 - gradient descent algorithm using Lagrangian in PCA-L1 paper
% X : PCA whitened data
% N : number of samples
% D : number of sources (sub-Gaussian only)

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
for ns=1:D
    w = randn(D,1);
    if ns > 1
        w = w - W*(W'*w);
    end
    w = w/norm(w);
    for iter = 1:max_iter
        G1 = zeros(D,1);
        for i=1:N
            x = X(:,i);
            a = w'*x;
            g1 = sign(a)*x; % for p=1
            G1 = G1+g1;
        end
        
        w_old = w;
        w = G1;
        % Gram-Schmidt orthogonalization
        if ns > 1
            w = w - W*(W'*w);
        end
        w = w/norm(w);
        % convergence check
        if norm(w-w_old) < Tol
            break;
        end
    end
    % concatenation of unmixing vector into the unmixing matrix
    W = [W,w];
    
end