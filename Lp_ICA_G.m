function W = Lp_ICA_G(X, D, N)
% Lp ICA - gradient descent algorithm
% X : PCA whitened data
% N : number of samples
% D : number of sources

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
for ns = 1:D
    
    w = randn(D,1);
    % orthogonalization by Gram-Schmidt
    if ns > 1
        w = w - W*(W'*w);
    end
    w = w/norm(w);
    
    % use log likelihood to find the optimal parameter p (super- or sub-G)
    p = 1;
    %p = plist(pind);
    
    % finding an unmixing vector
    for iter = 1:max_iter
        % Gradient; G1 = gradient when p=1, G3 = gradient when p=3
        G1 = zeros(D,1);
        G3 = zeros(D,1);
        
        % Log likelihood of 3 models
        LL_super = N*log(magnitude(1)); % super-G: p=1
        LL_gau = N*log(magnitude(2)); % G: p=2
        LL_sub = N*log(magnitude(3)); % sub-G: p=4
        
        for i=1:N
            x = X(:,i);
            a = w'*x;
            g1 = sign(a)*x; % for p=1
            g3 = sign(a)*power(abs(a),2)*x; % for p=3
            G1 = G1+g1;
            G3 = G3+g3;
            
            LL_super = LL_super - scale(1)*power(abs(a),1);
            LL_gau = LL_gau - scale(2)*power(abs(a),2);
            LL_sub = LL_sub - scale(3)*power(abs(a),3);
        end
        
        w_old = w;
        
        if LL_super > LL_sub
            w = w - alpha*G1;    % minimization for p=1
        else
            %w = w + alpha*G1;    % maximization for p=1
            %w = G1;      % maximization by Lagrange method (for p=1)
            w = w - alpha*G3;   % minimization for p=3
        end
        % Gram-Schmidt orthogonalization
        if ns > 1
            w = w - W*(W'*w);
        end
        w = w/norm(w);
        
        % convergence check
        if norm(w-w_old) < Tol
            break;
        end
    end  % end of iter
    %{
    % check whether the found source is super- or sub- G
    if (LL_super > LL_sub) & (LL_super > LL_gau)
        source_type(ns,:) = 'sup-G';
    elseif (LL_sub > LL_gau)
        source_type(ns,:) = 'sub-G';
    else
        source_type(ns,:) = 'Gauss';
    end
    %}
    % concatenation of unmixing vector into the unmixing matrix
    W = [W,w];
        
end % end of ns