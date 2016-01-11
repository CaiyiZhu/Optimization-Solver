function val=deri2(x)
%%% return the hessian matrix of the objective
load data/ABD
%%%A:100*1: the x coordinate of the 100 stores
%%%B:100*1: the y coordinate of the 100 stores
x1=x(1:50);
x2=x(51:100);
d1=zeros(50,1);%%%% the second order derivative of f over df/dx_idx_i, i=1,2,...,50
for i=1:50
    d1(i)=sum(D.*(x2(i)-B).^2.*((x1(i)-A).^2+(x2(i)-B).^2).^(-3/2));
end

d2=zeros(50,1); %%%% the second order derivative of f over df/dy_idy_i, i=1,2,...,50
for i=1:50
    d2(i)=sum(D.*(x1(i)-A).^2.*((x1(i)-A).^2+(x2(i)-B).^2).^(-3/2));
end

d3=zeros(50,1); %%%% the second order derivative of f over df/dx_i dy_i, i=1,2,...,50
for i=1:50
    d3(i)=-sum(D.*(x1(i)-A).*(x2(i)-B).*((x1(i)-A).^2+(x2(i)-B).^2).^(-3/2));
end

H=zeros(100,100);

for i=1:50
    H(i,i)=d1(i);
end
for i=1:50
    H(50+i,50+i)=d2(i);
end
for i=1:50
    H(i,50+i)=d3(i);
    H(50+i,i)=d3(i);
end

val=H;

end