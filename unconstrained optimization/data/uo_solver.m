function x_optimal = uo_solver(solver_method, data)

%% INPUTS      
%
% solver_method     string which contains the name of the four
%                   unconstrained optimization methods studied in class,
%                   specifying which method is to be used to solve the
%                   problem
%                       - 'steepest_descent'
%                       - 'newtons'
%                       - 'broyden'
%                       - 'alternate_broyden'
%
% data              data structure which contains all required parameters
%                   to solve unconstrained optimization problem
%                       - 'steepest_descent'
%                           (a) data.f = symbolic function to be optimized
%                           (b) data.x0 = initial estimate
%                           (c) data.e = stopping criteria
%                           (d) data.line_search = a string specifying
%                           line search method which should be used: 
%                               (1) 'dich': Use dichotomous search
%                                   data must then also contain:
%                                       data.a, data.b = search interval
%                                       data.e1, data.e2 = tolerances
%                               (2) 'golden': Use golden search
%                                   data must then also contain:
%                                       data.a, data.b = search interval
%                                       data.e2: tolerance
%                               (3) 'armijo': Use Armijo's rule
%                                   data must then also contain:
%                                       data.lambda0 = initial step size
%                                       data.alpha = Armijo's parameter
%                       - 'newtons'
%                           (a) data.f = symbolic function to be optimized
%                           (b) data.x0 = initial estimate
%                           (c) data.e = stopping criteria
%                       - 'broyden'
%                           (a) data.f = symbolic function to be optimized
%                           (b) data.x0 = initial estimate
%                           (c) data.e = stopping criteria
%                           (d) data.phi = scalar parameter
%                           (e) data.D0 = SPD matrix
%                           (f) data.line_search = a string specifying
%                           line search method which should be used: 
%                               (1) 'dich': Use dichotomous search
%                                   data must also contain:
%                                       data.a, data.b = search interval
%                                       data.e1, data.e2 = tolerances
%                               (2) 'golden': Use golden search
%                                   data must also contain:
%                                       data.a, data.b = search interval
%                                       data.e2: tolerance
%                               (3) 'armijo': Use Armijo's rule
%                                   data must also contain:
%                                       data.lambda0 = initial step size
%                                       data.alpha = Armijo's parameter
%                       - 'alternate_broyden'
%                           (a) data.f = symbolic function to be optimized
%                           (b) data.x0 = initial estimate
%                           (c) data.e = stopping criteria
%                           (d) data.phi = scalar parameter
%                           (e) data.B0 = SPD matrix
%                           (f) data.line_search = a string specifying
%                           line search method which should be used: 
%                               (1) 'dich': Use dichotomous search
%                                   data must also contain:
%                                       data.a, data.b = search interval
%                                       data.e1, data.e2 = tolerances
%                               (2) 'golden': Use golden search
%                                   data must also contain:
%                                       data.a, data.b = search interval
%                                       data.e2: tolerance
%                               (3) 'armijo': Use Armijo's rule
%                                   data must also contain:
%                                       data.lambda0 = initial step size
%                                       data.alpha = Armijo parameter
%                                       data.epsilon = Armijo parameter
%
%% MAIN FUNCTION
epsilon=data.e;
xk=data.x0;
line_search=data.line_search;
switch solver_method
    
    case 'steepest_descent' 
        re=[];
        % EXECUTE STEEPEST DESCENT ALGORITHM
        while (1)
            %%%% compute direction: gradient descent
            val=data.objective(xk)
            gradient=deri1(xk);
            d=-gradient/norm(gradient);
            re=[re;val];
            %%% stopping criteria
            if(norm(gradient)<epsilon*(1+abs(val)))
                break;
            end
        
            % For determining step size, call appropriate line search function
            switch line_search
                case 'armijo'
                    % Use Armijo's rule
                    lambdak = armijo(data.objective, data.deri1, xk, d, data.lambda0, data.alpha, data.e);
                case 'dich'
                    % Use Dichotomous search algorithm
                    lambdak = dich(data.objective, xk, d, data.a, data.b, data.e1, data.e2);
                case  'golden'
                    % Use golden search algorithm
                    lambdak = golden(data.objective, xk, d, data.a, data.b, data.e2);
                case  'fibonacci'
                    % Use fibonacci search algorithm
                    lambdak = fibonacci(data.objective, xk, d, data.a, data.b, data.e2);
                case  'bisection'
                    % Use bisection search algorithm
                    lambdak = bisection(data.deri1, xk, d, data.a, data.b, data.e2);
                case  'newton_search'
                    % Use newton's search algorithm
                    lambdak = newton_search(data.objective, data.deri1, data.deri2, xk, d, data.lambda0, data.e2);       
            end
        %%%% update x
            xk=xk+lambdak*d;
        end
        figure
        plot(1:length(re),re)
        
        
    case 'newtons'
        re=[]
        while (1)
            %%%% compute direction: gradient descent
            val=data.objective(xk)
            re=[re;val];
            gradient=deri1(xk);
            d=-inv(deri2(xk))*gradient;
            

            if(norm(gradient)<epsilon*(1+abs(val)))
                break;
            end
            lambdak=1;
            %%%% update x
            xk=xk+lambdak*d;
        end
        figure
        plot(1:length(re),re)

    case 'broyden'
        % EXECUTE BROYDEN FAMILY ALGORITHM
        D=data.D0;
        phi=data.phi;
        re=[];
        while (1)
            %%%% compute direction: gradient descent
            val=data.objective(xk)
            re=[re;val];
            gradient=deri1(xk);
            d=-D*deri1(xk);
            %%%
            if(norm(gradient)<epsilon*(1+abs(val)))
                break;
            end

            % For determining step size, call appropriate line search function
            switch line_search
                case 'armijo'
                    % Use Armijo's rule
                    lambdak = armijo(data.objective, data.deri1, xk, d, data.lambda0, data.alpha, data.e);
                case 'dich'
                    % Use Dichotomous search algorithm
                    lambdak = dich(data.objective, xk,d, data.a, data.b, data.e1, data.e2);
                case  'golden'
                    % Use golden search algorithm
                    lambdak = golden(data.objective, xk,d, data.a, data.b, data.e2);
                case  'fibonacci'
                    % Use fibonacci search algorithm
                    lambdak = fibonacci(data.objective, xk,d, data.a, data.b, data.e2);
                case  'bisection'
                    % Use bisection search algorithm
                    lambdak = bisection(data.deri1, xk,d, data.a, data.b, data.e2);
                case  'newton_search'
                    % Use newton's search algorithm
                    lambdak = newton_search(data.objective, data.deri1, data.deri2, xk, d, data.lambda0, data.e2);       
            end
            %%%% update x
            x_prev=xk;
            xk=xk+lambdak*d;
            %update D:
            p=xk-x_prev;
            q=deri1(xk)-deri1(x_prev);
            v=p/(p'*q)-D*q/(q'*D*q);
            C=p*p'/(p'*q)-D*q*q'*D/(q'*D*q)+ phi*q'*D*q*v*v';
            D=D+C;
        end
        figure
        plot(1:length(re),re)
        
    case 'alternate_broyden'
        
        % EXECUTE ALTERNATE BROYDEN FAMILY ALGORITHM
        phi=data.phi;
        B=data.B0;
        re=[];
        while (1)
            val=data.objective(xk)
            re=[re;val];
            gradient=deri1(xk);
            d=-inv(B)*deri1(xk);
            if(norm(gradient)<epsilon*(1+abs(val)))
                break;
            end
            % For determining step size, call appropriate line search function
            switch line_search
                case 'armijo'
                    % Use Armijo's rule
                    lambdak = armijo(data.objective, data.deri1, xk, d, data.lambda0, data.alpha, data.e);
                case 'dich'
                    % Use Dichotomous search algorithm
                    lambdak = dich(data.objective, xk,d, data.a, data.b, data.e1, data.e2);
                case  'golden'
                    % Use golden search algorithm
                    lambdak = golden(data.objective, xk,d, data.a, data.b, data.e2);
                case  'fibonacci'
                    % Use fibonacci search algorithm
                    lambdak = fibonacci(data.objective, xk,d, data.a, data.b, data.e2);
                case  'bisection'
                    % Use bisection search algorithm
                    lambdak = bisection(data.deri1, xk,d, data.a, data.b, data.e2);
                case  'newton_search'
                    % Use newton's search algorithm
                    lambdak = newton_search(data.objective, data.deri1, data.deri2, xk, d, data.lambda0, data.e2);       
            end
            %%%% update x
            x_prev=xk;
            xk=xk+lambdak*d;
            %update D:
            p=xk-x_prev;
            q=deri1(xk)-deri1(x_prev);
            v=q/(q'*p)-B*p/(p'*B*p);
            C=q*q'/(q'*p)-B*p*p'*B/(p'*B*p)+ phi*q'*B*q*v*v';
            B=B+C;
        end
        figure
        plot(1:length(re),re);
        
    otherwise
        warning('Solver method not specified correctly');
end
x_optimal=xk;

end

%% LINE SEARCH ALGORITHMS

function lambdak = armijo(objective, deri1, xk, d, lambda0, alpha, epsilon)
% Armijo line search
gradient=deri1(xk);
t=lambda0; % step size
while 1
    if objective(xk+t*d) <= (objective(xk)+ t*epsilon*gradient'*d)
        break
    else
        t=t/alpha;
        if t<10e-20
            break
        end

    end
end
lambdak=t;
end

function lambdak = dich(objective, xk, d, a, b, e1, e2)
% Dichotomous search algorithm
%d is the update direction
while 1
    if (b-a)<2*e2
        break
    else
        lambda1=(a+b)/2-e1/2;
        lambda2=(a+b)/2+e1/2;
        if objective(xk+lambda1*d)>objective(xk+lambda2*d)
            a=lambda1;
        else
            b=lambda2;
        end
    end
end

lambdak=(a+b)/2;

end


function lambdak = golden(objective, xk, d, a, b, e2)
%golden search algorithm
%d is the update direction
while 1
    if (b-a)<e2
        break
    else
        lambda1=a+0.382*(b-a);
        lambda2=a+0.618*(b-a);
        if objective(xk+lambda1*d)>objective(xk+lambda2*d)
            a=lambda1;
        else
            b=lambda2;
        end
    end
    a;
    b;
end

lambdak=(a+b)/2;

end



function lambdak = fibonacci(objective, xk, d, a, b, e2)
%Fibonacci search
%d is the update direction

[n F]=find_n_FS((b-a)/e2);
% n is the number of updates
% F is a vector store the value of fabonacci series[f0,f1,...,fn]
p = a + (F(n-2)/F(n))*(b-a);
q = a + (F(n-1)/F(n))*(b-a);

fp = objective(xk+p*d);
fq = objective(xk+q*d);
% Step 3: perform repeated reductions of search interval using k as an
% iteration index

for k = n-1:-1:3
    if (fp <= fq)
        b = q;
        q = p;
        fq = fp;    % given that we set q = p, we also need to set fq = fp
                    % this is one of the reasons why the Fibonacci search is "efficient": 
                    % it re-uses the result of an f-calculation from the
                    % previous iteration.
        p = b - (F(k-1)/F(k))*(b-a);  
        fp = objective(xk+p*d);
    else    
        a = p;
        p = q;
        fp = fq;    % given that we set p = q, we also need to set fp = fq
        
        q = a + (F(k-1)/F(k))*(b-a);
        fq = objective(xk+q*d);
        
    end
end
lambdak=(a+b)/2;
end

function lambdak = bisection(deri1, xk, d, a, b, e2)
%bisection search, the first derivative of theta(lambda) over lambda is
%necessary
%d is the update direction
while 1
    if (b-a)<e2
        break
    else
        lambda=(a+b)/2;
        if d'*deri1(xk+lambda*d)>0 % if the increment is positive, then decrease the stepsize
            b=lambda;
        else
            a=lambda;
        end
    end
end
lambdak=(a+b)/2;
end


function lambdak = newton_search(objective, deri1, deri2, xk, d, lambda0, e2)
% Newton's line search algorithm
%d is the update direction
p=0;
q=lambda0;
while 1
    if abs(d'*deri1(xk+q*d))<e2*(1+abs(objective(xk+q*d)))%same stopping criteria as we used in newton's method
        break
    else
        % theta_1st_gradient=d'*grad(f(x+lambda*d)),
        % theta_2nd=d'*H(f(x+lambda*d))*d
        p=q;
        q=p-(d'*deri1(xk+p*d))/(d'*deri2(xk+p*d)*d);
    end
end
lambdak=q;
end