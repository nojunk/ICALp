function X = whitening( X0 )
%UNTITLED2 이 함수의 요약 설명 위치
%   자세한 설명 위치

[d,n] = size(X0);
[U,L,V] = svd(X0);
sv = diag(L);
invL = diag(1./sv)*sqrt(n);

X = invL*U'*X0;
end

