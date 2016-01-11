function val=deri1(x)
%%% return the first order derivative of the objective
load data/ABD

val=zeros(100,1);%% the first derivative of the objective over x and y, the first 50 entries are over x and the next entries are for y
x1=x(1:50);
x2=x(51:100);
for i=1:50
    for j=1:100
        val(i)=val(i)+D(j)*(x1(i)-A(j))/sqrt((x1(i)-A(j))^2+(x2(i)-B(j))^2);
    end
end

for i=1:50
    for j=1:100
        val(i+50)=val(i+50)+D(j)*(x2(i)-B(j))/sqrt((x1(i)-A(j))^2+(x2(i)-B(j))^2);
    end
end
end
