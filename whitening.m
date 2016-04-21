function X = whitening( X0 )
%UNTITLED2 �� �Լ��� ��� ���� ��ġ
%   �ڼ��� ���� ��ġ

[d,n] = size(X0);
[U,L,V] = svd(X0);
sv = diag(L);
invL = diag(1./sv)*sqrt(n);

X = invL*U'*X0;
end

