function [ P ] = pol_n_permute( B )
%UNTITLED3 이 함수의 요약 설명 위치
%   자세한 설명 위치
[s,d] = size(B);

T = abs(B);    % temporary
P = zeros(s,d);

for i=1:s
    [val, r_ind] = max(T);
    [val, c(i)] = max(val);
    r(i) = r_ind(c(i));
    
    T(r(i),:) = 0;
    T(:,c(i)) = 0;
    
    P(r(i),c(i)) = sign(B(r(i),c(i)));
end

end

