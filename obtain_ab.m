function [alpha, beta] = obtain_ab(p_list, load_option)
%obtain_ab: obtain alpha and beta for p(x) = alpha* exp(-beta*|x|^p)
%   Nojun Kwak (nojunk@snu.ac.kr)
%   Sep. 11, 2015

dx = 0.001;

if nargin < 2 
   load_option = false;
end

if load_option == true
    load('alpha_beta.mat');

else
for q=1:length(p_list)
    
    b(q) = 0.0;
    for i=-10:dx:10
        f = exp(-power(abs(i),p_list(q)));
        b(q) = b(q) + f*dx;
    end
    
    var(q) = 0.0;
    for i=-10:dx:10
        f = exp(-power(abs(i),p_list(q)))/b(q);
        var(q) = var(q) + power(i,2)*f*dx;
    end
    
    %p = exp(-power(abs(i),q))/b;
    %plot(i,p)
    
    s(q) = sqrt(var(q));
    alpha(q) = power(s(q),p_list(q));
    beta(q) = s(q)/b(q); 
    
%     p = exp(-alpha(q)*power(abs(x),q_list(q)))*beta(q);
%     data = [data;p];
%     hold on;
%     %plot(x,p,'k');
    
%     c(q) = 0;
%     for i=-10:dx:10
%         p = exp(-alpha(q)*power(abs(i),q_list(q)))*beta(q);
%         c(q) = c(q) + p*dx;
%     end
%     
end

end

end
    
