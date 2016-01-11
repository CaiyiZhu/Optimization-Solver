clear

% Define the objective function

% %% large problem
syms x;
data.n = 84;
a = [24,20,6,33,50,34,22,33,30,38,33,26];
b = [49,15,17,50,40,40,32,41,23,42,7,26];
d = [27,45,33,24,5,49,32,8,20,45,27,28];
c = [41,27,30,12,6,33];
na=length(a); % na=12
nc=length(c); % nc=6
data.f = @(x) 0;
for i = 1:nc
    for j = 1:na
        data.f = @(x) data.f(x) + x( na + (i-1)*na + j).*sqrt((x(i) - a(j)).^2 + ( x(nc+i) - b(j)).^2);
        data.f = @(x) data.f(x) + x( na + (i-1)*na + j).*((x(i) - a(j)).^2 + (x(6+i) - b(j)).^2);
    end
end
       
% Define the inequality constraints, gi, i = 1,...,m
for i = 1:nc
    data.g{i} = @(x) 0;
    for j = 1:na
        data.g{i} = @(x) data.g{i}(x) + x(na + (i-1)*na + j);
    end
    data.g{i} = @(x) data.g{i}(x) - c(i);
end
for i=na+1:na+72
    data.g{i-nc} =@(x) - x(i);
end
m= nc + 72; % number inequality constraints

% Define the equality constraints, hi, i = 1,...l
for j = 1:na
    data.h{j} = @(x) 0;
    for i = 1:nc
        data.h{j} = @(x) data.h{j}(x) + x(na + (i-1)*na + j);
    end
    data.h{j} = @(x) data.h{j}(x) - d(j);
end
l=nc; % number of equality constraints



%% LARGE PROBLEM, The Other Student
% Comment out if running small problem

% % Define the objective function
% syms x;
% data.n = 30;
% data.f = @(x) 0;
% b = 0.1:0.1:data.n/10;
% for i = 1:data.n
%     data.f = @(x) data.f(x) + x(i).^2 + b(i).*(10 - x(i)).^4;
% end
% 
% % Define the inequality constraints, gi, i = 1,...,m
% for i = 1:data.n
%     data.g{i} = @(x) -x(i)^2 + 4;
% end
% for i = 1:data.n
%     data.g{data.n + i} = @(x) -(10 - x(i)).^2 + 4;
% end
% m=2*data.n;
% 
% % Define the equality constraints, hi, i = 1,...l
% data.h = {};
% 
% l=0;




% %% small problem
% data.f = @(x) (x(1)-2)^2 + (x(2)-2)^2;
% data.n=2;
% data.g={};
% m=1; % number of inequality constraints
% for i =1:m
%     data.g{i}= @(x) x(1)^2/4 + x(2)^2 -1;
% end
% 
% 
% data.h={};
% l=1; % number of equality constraints
% for i=1:l
%     data.h{i}=@(x) x(1)-2*x(2)+1;
% end


%% initialization 
% options = optimset('Display','iter','TolX',1e-20);
% optnew = optimset(options,'TolX',1e-10);
data.c0=100;
data.beta=5; % for penalty method, beta is greater than 1, and for barrier method, beta should be less than 1
% data.x0=rand(data.n,1);
data.epsilon=1e-10;
% data.sym=x;
data.u0=rand(m,1);
data.v0=rand(l,1);
data.VIOLx0=10000;
data.c1=1;
data.B0=eye(data.n);
data.m=m;
data.l=l;
% 	                    - 'penalty'
%                       - 'barrier', need to transform the equality
%                       constraint
%                       - 'ALAG' 
%                       - 'SQP'
%                       - 'SQPBFGS'
%                       - 'MSQP'

% % penalty method and ALAG method
% data.c0=10;
% data.x0=ones(data.n,1); % for penalty method, the initial point should be infeasible, for barrier method, the initial point should be feasible
% data.beta=5; % for penalty method, beta is greater than 1, and for barrier method, beta should be less than 1
% solver_method='penalty';   % or 'ALAG', penalty method works
% 
% % barrier method, need to convert the equality constraint into inquality
% %constant
% data.c0=10;
% data.x0=0.1*ones(data.n,1); % for barrier method, the initial point should be feasible 
% data.beta=0.5; % for penalty method, beta is greater than 1, and for barrier method, beta should be less than 1
% solver_method='barrier';   %%% the result returned by barrier method is the real optimal solution

%SQP, SQPBFGS, MSQP
data.x0=rand(data.n,1); % close initial data for small problem, to gen
% data.x0=10*rand(data.n,1); % far initial data for small problem
data.x0=20*ones(data.n,1);
solver_method='SQP';  

[x_optimal, obj_value, runtime] = co_solver(solver_method, data)



