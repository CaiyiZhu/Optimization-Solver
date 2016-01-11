function val=deri1(x)
% f(x1,x2)=(x1+x2)^2+exp^(x1+x2)+sin(x1+x2), return the first order of
% derivative
% s=sum(x);
% val=(2*(s)+exp(s)+cos(s))*ones(length(x),1);

%update: f(x1,x2)=(x1-1)^2+2*(x2-2)^2 + sin(x1+x2)
x1=x(1);x2=x(2);
val=[2*(x1-1)+cos(x1+x2);2*(x2-2)+cos(x1+x2)];