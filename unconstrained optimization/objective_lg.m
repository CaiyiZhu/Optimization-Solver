function val=objective(x)
%x:100*1: the first 50 entries of x is x coordinates of the 50 warehouse, the last 50 entries of x is
%the y coordinates of the 50 warehouse
load data/ABD

x1=x(1:50);
x2=x(51:100);

val=0;
%c=1; % constant in this problem
for i=1:50
    for j=1:100
        val=val+D(j)*sqrt((x1(i)-A(j))^2+(x2(i)-B(j))^2);
    end
end
%val=val*c;

end
