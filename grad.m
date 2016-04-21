function [fpval] = grad(fval, dp)
%grad: taking the gradient (derivative) of the function value list fval
%   Nojun Kwak (nojunk@snu.ac.kr)
%   Sep. 11, 2015

n = length(fval);
fpval(1) = (fval(2) - fval(1))/dp;
fpval(n) = (fval(n) - fval(n-1))/dp;
for i=2:n-1
    fpval(i) = (fval(i+1) - fval(i-1))/2/dp;
end

end
