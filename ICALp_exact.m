function [ W, Y, P, LL ] = ICALp_exact( D, s, init_p, whitened )
%ICALp_exact Independent component analysis by Lp-norm minimization 
%   Nojun Kwak (nojunk@snu.ac.kr)
%   Sep. 09, 2015
%
%   D: d * n data matrix, assumed to be whitened
%   s: no. of independent components
%   init_p: initial value of p, default = 2.0
%   W: d * s unmixing matrix 
%   Y: s * n independents sources Y = W'*X
%   whitened: whether X is whitened or not (default: whitened = true)
%   P: optimal p (in Lp norm)
%   LL: log-likelihood


X = D;
[d,n] = size(X);

if nargin < 2
    s = d;       
end

if nargin < 3
    init_p = 2.0;   
end

if nargin < 4
    whitened = true;
end

if whitened == false
    [X, Z] = whitening(D);
end



% obtain alpha and beta for various values of p
stepp = 0.01;
p_list = 0.01:stepp:10;
[alpha, beta] = obtain_ab(p_list,true);
dalpha = grad(alpha, stepp);
dbeta = grad(beta,stepp);


MAX_ITER = 1000;
LL = zeros(s,MAX_ITER);
W = zeros(d,s);
for i=1:s
    % initialization 
    p = init_p;    % start from Gaussian
    w = randn(d,1);   
    if i > 1
        X = X - W(:,i-1)*(W(:,i-1)'*X);
        w = w - W(:,1:i-1)*(W(:,1:i-1)'*w);    
    end
    if i==s
        converged = true;
    else
        converged = false;
    end
    
    w = w/norm(w);
      
    pstep = 0.0;    
    wstep = zeros(d,1);    
    rho_base = 0.1;   
    
    LL_past = -1E10;
    n_iter = 0;    
    while converged == false && n_iter < MAX_ITER      
        rho = rho_base*(p^2);
        n_iter = n_iter+1;        
        pind = round(p/stepp);
        y = X'*w;
        absy = abs(y);
        sgny = sign(y);
        labsy = log(absy);
        pabsy = power(absy,p);
        s1 = sum(pabsy);
        s2 = sum(labsy.*pabsy);
        
        if (p > 0.7) && (p < 5)            
            dp = dbeta(pind)/beta(pind) - dalpha(pind)*s1/n - alpha(pind)*s2/n;                
            p = p + rho*dp;    
        else
            dp = 0;
        end
        
        s3 = zeros(d,1);
        for j=1:n
            s3 = s3 + power(absy(j),p-1) * sgny(j) * X(:,j);
        end        
        dw = -alpha(pind)*s3/n;
        
        LL(i,n_iter) = n*log(beta(pind)) - alpha(pind)*s1;        
        
        w_old = w;
        w = w + rho*dw;
        if i > 1            
            w = w - W(:,1:i-1)*W(:,1:i-1)'*w;
        end
        w = w/norm(w); 
        
        if (norm([w_old-w;dp]) < 1E-10)
            converged = true;                    
        end
        
    end
    W(:,i) = w;    
    P(i) = p;
end
Y = W'*D;

end

