function val=deri2(x)
% f(x1,x2)=(x1+x2)^2+exp^(x1+x2)+sin(x1+x2), return the 2nd order of
% derivative, i.e., the hessian matrix
% s=sum(x);
% val=(2+exp(s)-sin(s))*ones(length(x),length(x));

%update: f(x1,x2)=(x1-1)^2+2*(x2-2)^2 + sin(x1+x2)
x1=x(1);x2=x(2);
val=[2 -sin(x1+x2);-sin(x1+x2) 2];