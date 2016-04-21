function [W2, W1]=FICA1(x, ini, g)
%algorithm 1FICA
% 
% version: 1.0  release: 19.1.2007
%
% copyright: Zbynek Koldovsky, Petr Tichavsky



epsilon=0.000001;
[dim N]=size(x);
MaxIt=100;

repeat=1;
rot2d=[1/sqrt(2) 1/sqrt(2);-1/sqrt(2) 1/sqrt(2)];
SaddleTest=true;

if nargin<3
    g='tanh';
end

if nargin<2
   ini=randn(dim,dim);
end

test_of_saddle_points_nonln='tanh';



for j=1:dim
  x(j,:)=x(j,:)-mean(x(j,:));
end 



%preprocessing
C = cov(x');


%whitened data
x = C^(-1/2)*x;


%%% FastICA Symmetric approach with the test of saddle points
W=ini;
W=W*real(inv(W'*W)^(1/2));
while repeat
crit=0; NumIt=0;
 while (1-min(crit)>epsilon && NumIt<MaxIt)% && sum(double(changed>10))<2)
   Wold=W;
   switch g
    case {'tanh','biga'}
     hypTan = tanh(x'*W);
     W=x*hypTan/N-ones(dim,1)*sum(1-hypTan.^2).*W/N;
    case 'pow3'
     W=(x*(pwr(x'*W,3)))/N-3*W;
    case 'gaus'
     U=x'*W;
     Usquared=U.^2;
     ex=exp(-Usquared/2);
     gauss=U.*ex;
     dGauss=(1-Usquared).*ex;
     W=x*gauss/N-ones(dim,1)*sum(dGauss).*W/N;
   end
   W=W*real(inv(W'*W)^(1/2));
   crit=abs(sum(W.*Wold));
   NumIt=NumIt+1;
 end %while iteration
 repeat=0;
%%%The test of saddle points
 if SaddleTest
  SaddleTest=false; %%The test could be done only one times
  u=x'*W;
  switch test_of_saddle_points_nonln
        case 'tanh'
            table1=(mean(log(cosh(u)))-0.37456).^2;
        case 'gaus'
            table1=(mean(ex)-1/sqrt(2)).^2;
        case 'pow3'
            table1=(mean((pwr(u,4)))-3).^2;
  end
  rotated=zeros(1,dim);
  checked=1:dim; 
  for i=checked
    for j=checked(checked>i)
        if (~rotated(i) && ~rotated(j))
        h=[u(:,i) u(:,j)]*rot2d;
        switch test_of_saddle_points_nonln
            case 'tanh'
                ctrl=(mean(log(cosh(h)))-0.37456).^2;
            case 'gaus'
                ctrl=(mean(exp(-h.^2/2)-1/sqrt(2))).^2;
            case 'pow3'
                ctrl=(mean((pwr(h,4)))-3).^2;
        end
        if sum(ctrl)>table1(i)+table1(j)
            %bad extrem indicated
            rotated([i j])=1; %do not test the rotated signals anymore
            W(:,[i j])=W(:,[i j])*rot2d;
            repeat=1; %continue in iterating - the test of saddle points is positive
            MaxIt=30;
        end
        end
    end
  end
 end %if SaddleTest
end %while repeat


 
W1=W'*C^(-1/2);


%one-unit stage (for all sources)
for k=1:dim
  w=W(:,k);
  w=w/norm(w);
  wold=zeros(dim,1);
  NumIt=0;
 while (abs(w'*wold)<1-epsilon && NumIt<100 && ...
        (abs((W(:,k)/norm(W(:,k)))'*(w/norm(w)))>0.90))
   wold=w;
   u=w'*x;
    switch g
     case 'tanh'
      hyptan=tanh(u);
      w=(x*hyptan' - sum(1-hyptan.^2)'*w)/N;
     case 'gaus'
      gauss=u.*exp(-u.^2/2);
      dGauss=(1-u.^2).*exp(-u.^2/2);
      w=(x*gauss'-sum(dGauss)*w)/N;
     case 'biga'
      u=u';   
      m=mean(abs(u)); %estimate centers of distribution's 
      e=sqrt(1-m^2); %then their variance is..., because the overall variance is 1
      if e<=0.05, e=0.05; m=sqrt(1-e^2); end %due to stability
      uplus=u+m; uminus=u-m;
      expplus=exp(-uplus.^2/2/e^2);
      expminus=exp(-uminus.^2/2/e^2);
      expb=exp(-(u.^2+m^2)/e^2);
      gg=-(uminus.*expminus + uplus.*expplus)./(expplus+expminus)/e^2;
      gprime=-(e^2*(expplus.^2+expminus.^2)+(2*e^2-4*m^2)*expb)./(expplus+expminus).^2/e^4;
      w=x*gg/N-mean(gprime)*w;                   
     case 'pow3'
      w=(x*(pwr(u,3))')/N-3*w;
    end
    w=w/norm(w);
    NumIt=NumIt+1;
 end
 if (abs((W(:,k)/norm(W(:,k)))'*(w/norm(w)))>0.95)
     W(:,k)=w;
 end
end

W2=W'*C^(-1/2);




function x=pwr(a,n)
x=a;
for i=2:n
    x=x.*a;
end
