function val=objective(x)
% f(x1,x2)=(x1+x2)^2+exp^(x1+x2)+sin(x1+x2), give an input calcalate the value of the function

% s=sum(x);
% val=(s)^2+exp(s)+sin(s);

% update: f(x1,x2)=(x1-1)^2+2*(x2-2)^2 +sin(x1+x2)
x1=x(1);x2=x(2);
val=(x1-1)^2+(x2-2)^2 +sin(x1+x2) ;%;+x1*x2;