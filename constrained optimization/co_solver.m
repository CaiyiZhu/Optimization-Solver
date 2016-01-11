function [x_optimal, obj_value, runtime] = co_solver(solver_method, data)
%% INPUTS      
%
% solver_method     string which contains the name of the six constrained 
%                   optimization methods studied in class, specifying which
%                   method is to be used to solve the problem
%                       - 'penalty'
%                       - 'barrier'
%                       - 'ALAG'
%                       - 'SQP'
%                       - 'SQPBFGS'
%                       - 'MSQP'
%
% data              data structure which contains all required parameters
%                   to solve unconstrained optimization problem
%                   - In all cases, data must contain the following:
%                       - data.f = symbolic objective function
%                       - data.sym = vector of symbols used to create f
%                       - data.g = vector containing inequality constraints
%                       - data.h = vector containing equality constraints
%                       - data.epsilon = termination criteria
%                   - For each of the respective solver methods, the
%                   following must also be included
%                       - 'penalty'
%                           (a) data.c0 = initial penalty parameter
%                           (b) data.beta = acceleration parameter
%                       - 'barrier'
%                           (a) data.c0 = initial penalty parameter
%                           (b) data.beta = acceleration parameter
%                       - 'ALAG'
%                           (a) data.x0 = initial solution
%                           (b) data.u0 = inequality constraint Lagrangian
%                           (c) data.v0 = equality constraint Lagrangian
%                           (d) data.c1 = initial penalty parameter
%                           (e) data.beta = acceleration parameter
%                           (f) data.VIOLx0 = violations at x0
%                       - 'SQP'
%                           (a) data.x0 = initial solution
%                           (b) data.u0 = inequality constraint Lagrangian
%                           (c) data.v0 = equality constraint Lagrangian
%                       - 'SQPBFGS'
%                           (a) data.x0 = initial solution
%                           (b) data.u0 = inequality constraint Lagrangian
%                           (c) data.v0 = equality constraint Lagrangian
%                           (d) data.B0 = initial Hessian approximation
%                       - 'MSQP'
%                           (a) data.x0 = initial solution
%                           (b) data.u0 = inequality constraint Lagrangian
%                           (c) data.v0 = equality constraint Lagrangian
%                           (d) data.B0 = initial Hessian approximation
%                           (e) data.c = convergence parameter
%
%% NOTES
%
% If gradients and/or Hessian matrices are required for algorithms, if you
% want to save time, you can use the jacobian() and/or hessian() functions
% to generate the symbolic representations at the beginning of the 
% algorithm (i.e. after case 'solver_method'). Then, use matlabFunction() 
% to convert these representations into a function. For example:
%       syms x
%       f = @(x) x^2;  % wrong message, no @(x), f should be a symbolic
%       function, not an anounymous function, and is only useful when x is
%       a scalar
%       grad_sym = jacobian(f);
%       grad_f = matlabFunction(grad_sym);
% This produces a MATLAB function grad_f = @(x) 2*x, as desired. So typing
% grad_f(2) will produce the answer 4. You can then use grad_f as a
% function throughout the algorithm. This MUST be done in co_solver, and
% not in run_file and passed as an element in data, as calculating the
% gradient and Hessian once is part of the algorithm.
%
% To solve the unconstrained subproblems in penalty, barrier, and ALAG,
% please use MATLAB's unconstrained solver, fminunc, details of which can
% be found at http://www.mathworks.com/help/optim/ug/fminunc.html.
%
% To solve the QP problem in SQP, SQPBFGS, and MSQP, please use MATLAB's QP
% solver, quadprog, details of which can be found at 
% http://www.mathworks.com/help/optim/ug/quadprog.html.
%% MAIN FUNCTION
% syms x
% f=@(x) (x(1)-2)^2 + (x(2)-2)^2;
% h=@(x) x(1)^2/4 + x(2)^2 -1;
% g=@(x) x(1) - 2*x(2);
tic;
fvalAll=[];
switch solver_method % 
    
    case 'penalty'
        % PENALTY ALGORITHM HERE
        EPSILON=data.epsilon;
        c=data.c0;
        beta=data.beta;
        q=2;
        %x=data.sym;
        p=@(x) 0;
        for i=1:data.m
            p =@(x) p(x) + max(0,data.g{i}(x))^q;
        end
        for i=1:data.l
            p =@(x) p(x) + abs(data.h{i}(x))^q;
        end
         
        
        k=1;
        x=data.x0;
        while c < 1e10
            L=@(x) data.f(x) + c*p(x);
            c
            [x,fval]=fminunc(L,x);
            fvalAll=[fvalAll data.f(x)];
            if p(x) < EPSILON
                break
            end
            c = beta*c;
            k = k+1;
%             penalty=p(x)
%             fval
            pause(1)
        end
%         p(x)
    case 'barrier'
        % BARRIER ALGORITHM HERE
        EPSILON=data.epsilon;
        c=data.c0;
        beta=data.beta;
        b=@(x) 0;
        for i=1:data.m
            b =@(x) b(x) - log(-data.g{i}(x));
        end
%         b = @(x) max(0,b(x)); % the penalty can not be negative
        x=data.x0;
        k=0;
        
        while c > 1e-20
            L =@(x) data.f(x) + c*b(x);
            [x,fval]=fminunc(L,x);
            fvalAll=[fvalAll data.f(x)]
            %if b(x) < EPSILON
               % break
            %else
                c = beta*c;
                k = k+1;
            %end
        end
        b(x)
    case 'ALAG'
        % ALAG ALGORITHM HERE

        EPSILON=data.epsilon;
        volPrev=data.VIOLx0;
        
        c=data.c1;
        u=data.u0;
        v=data.v0;
        beta=data.beta;

        
        
        x=data.x0;
        k=0;
        while c < 1e10
            c
            % define the function
            ALAG=@(x) data.f(x);
            for i=1:data.l
                ALAG =@(x) ALAG(x) + v(i)*data.h{i}(x) + c*(data.h{i}(x))^2;
            end
            for i=1:data.m
                ALAG =@(x) ALAG(x) + c*(max(0,data.g{i}(x)+u(i)/(2*c))^2) - u(i)^2/(4*c);
            end
            
            [x,fval]=fminunc(ALAG,x);
            fvalAll=[fvalAll data.f(x)]
            vol=0;
            for i=1:data.m
                if data.g{i}(x) > vol
                    vol=data.g{i}(x);
                end
            end
            for i=1:data.l
                if abs(data.h{i}(x)) > vol
                    vol=abs(data.h{i}(x));
                end
            end
            if vol < EPSILON
                break % violation is samll, terminate the code
            elseif vol < 0.25*volPrev
                vol
                volPrev
                % outer loop, update u and v
                for i=1:data.m
                    u(i) = max(0,u(i)+2*c*data.g{i}(x));
                end
                for i=1:data.l
                    v(i) = v(i) + 2*c*data.h{i}(x);
                end
%                 volPrev=data.VIOLx0;
                volPrev = vol;
%                 c=data.c1;
                u=data.u0;
                v=data.v0;
            else
                % continue inner loop, beta should be greater than 0
                 c = c*beta;
                 volPrev=vol;
            end
            k = k+1;
        end
                    
    case 'SQP'
        % SQP ALGORITHM HERE, the algorithm works for the small problem
        EPSILON=data.epsilon;
        x=data.x0;
        u=data.u0;
        v=data.v0;
        n=length(x);
        m=length(u);
        l=length(v);
        k=0;
        n
        %%% 
        sym_list=[];
        for  i=1:n
            sym_list=[sym_list sym(sprintf('x%d',i))];
        end
        % grad and hessian of f
        f_grad = gradient(data.f(sym_list),sym_list);
        f_hess = hessian(data.f(sym_list),sym_list);

        % grad and hessian of g
        for i=1:m
            g_grad{i} = gradient(data.g{i}(sym_list),sym_list);
            g_hess{i} = hessian(data.g{i}(sym_list),sym_list);
        end
        % grad and hessian of h
        for i=1:l
            h_grad{i} = gradient(data.h{i}(sym_list),sym_list);
            h_hess{i} = hessian(data.h{i}(sym_list),sym_list);            
        end

        % to evaluate the gradient and hessian,
        xk=x;
        while k<1000
            k
            xk
            f_grad_eval=eval(sym(subs(f_grad,sym_list,sym(xk))));
            f_hess_eval=eval(sym(subs(f_hess,sym_list,sym(xk))));
            for i=1:m
                g_grad_eval{i}=eval(sym(subs(g_grad{i},sym_list,sym(xk))));
                g_hess_eval{i}=eval(sym(subs(g_hess{i},sym_list,sym(xk))));
            end
            for i=1:l
                h_grad_eval{i}=eval(sym(subs(h_grad{i},sym_list,sym(xk))));
                h_hess_eval{i}=eval(sym(subs(h_hess{i},sym_list,sym(xk))));
            end
        
            % calculate the H_quad, f_quad
            H_quad=f_hess_eval;
            for i=1:m
                H_quad = H_quad + u(i)*g_hess_eval{i};
            end
            for i=1:l
                H_quad = H_quad + v(i)*h_hess_eval{i};
            end
            f_quad=f_grad_eval;
            
            A_quad=[];
            b_quad=[];
            for i=1:m
                A_quad = [A_quad g_grad_eval{i}];
                b_quad = [b_quad; -data.g{i}(xk)];
            end

            
            Aeq_quad=[];
            beq_quad=[];
            for i=1:l
                Aeq_quad = [Aeq_quad h_grad_eval{i}];
                beq_quad = [beq_quad; -data.h{i}(xk)];
            end
            
            [d,fval,exitflag,output,lambda]=quadprog(H_quad,f_quad,A_quad',b_quad,Aeq_quad',beq_quad);
            d
            xk = xk + d;
            fvalAll=[fvalAll data.f(xk)]
            u = lambda.ineqlin;
            v = lambda.eqlin;
            
            if norm(d) < EPSILON
                break
            end
            
        end
        x=xk;  %%% yes, it works    






        
    case 'SQPBFGS'
        % SQPBFGS ALGORITHM HERE, the algorithm works for the small problem
        EPSILON=data.epsilon;
        x=data.x0;
        u=data.u0;
        v=data.v0;
        n=length(x);
        m=length(u);
        l=length(v);
        B=data.B0;
        k=0;
        
        %%% 
        sym_list=[];
        for  i=1:n
            sym_list=[sym_list sym(sprintf('x%d',i))];
        end
        % grad and hessian of f
        f_grad = gradient(data.f(sym_list),sym_list);
        f_hess = hessian(data.f(sym_list),sym_list);

        % grad and hessian of g
        for i=1:m
            g_grad{i} = gradient(data.g{i}(sym_list),sym_list);
            g_hess{i} = hessian(data.g{i}(sym_list),sym_list);
        end
        % grad and hessian of h
        for i=1:l
            h_grad{i} = gradient(data.h{i}(sym_list),sym_list);
            h_hess{i} = hessian(data.h{i}(sym_list),sym_list);            
        end
        
%         % grad and hessian of L
%         L=@(x) data.f(x);
%         for i=1:m
%             L = @(x) L(x) + u(i)*data.g{i}(x);
%         end
%         for i=1:l
%             L = @(x) L(x) + v(i)*data.h{i}(x);
%         end
%         L_grad=gradient(L(sym_list),sym_list);
%         L_hess=hessian(L(sym_list),sym_list);

        % to evaluate the gradient and hessian,
        xk=x;
        xk_prev=xk;
        while k<1000
            k = k+1
            f_grad_eval=eval(sym(subs(f_grad,sym_list,sym(xk))));
            f_hess_eval=eval(sym(subs(f_hess,sym_list,sym(xk))));
            for i=1:m
                g_grad_eval{i}=eval(sym(subs(g_grad{i},sym_list,sym(xk))));
                g_hess_eval{i}=eval(sym(subs(g_hess{i},sym_list,sym(xk))));
            end
            for i=1:l
                h_grad_eval{i}=eval(sym(subs(h_grad{i},sym_list,sym(xk))));
                h_hess_eval{i}=eval(sym(subs(h_hess{i},sym_list,sym(xk))));
            end
        
            % calculate the H_quad, f_quad
            H_quad=B;
            f_quad=f_grad_eval;
            
            A_quad=[];
            b_quad=[];
            for i=1:m
                A_quad = [A_quad g_grad_eval{i}];
                b_quad = [b_quad; -data.g{i}(xk)];
            end

            
            Aeq_quad=[];
            beq_quad=[];
            for i=1:l
                Aeq_quad = [Aeq_quad h_grad_eval{i}];
                beq_quad = [beq_quad; -data.h{i}(xk)];
            end
            
            %update the hessian matrix
            L_grad_quad = f_grad_eval;
            for i=1:m
                L_grad_quad = L_grad_quad + u(i)*g_grad_eval{i};
            end
            for i=1:l
                L_grad_quad = L_grad_quad + v(i)*h_grad_eval{i};
            end
%             L_grad_eval_prev=eval(sym(subs(L,sym_list,sym(xk))));
            
            [d,fval,exitflag,output,lambda]=quadprog(H_quad,f_quad,A_quad',b_quad,Aeq_quad',beq_quad);
            d
            
            xk_prev=xk;
            xk = xk + d;
            fvalAll=[fvalAll data.f(xk)]
            u = lambda.ineqlin;
            v = lambda.eqlin;
            
            if norm(d) < EPSILON
                break
            end
            
            
            f_grad_eval_new=eval(sym(subs(f_grad,sym_list,sym(xk))));
            for i=1:m
                g_grad_eval_new{i}=eval(sym(subs(g_grad{i},sym_list,sym(xk))));
            end
            for i=1:l
                h_grad_eval_new{i}=eval(sym(subs(h_grad{i},sym_list,sym(xk))));
            end
            
            L_grad_quad_new = f_grad_eval_new;
            for i=1:m
                L_grad_quad_new = L_grad_quad_new + u(i)*g_grad_eval_new{i};
            end
            for i=1:l
                L_grad_quad_new = L_grad_quad_new + v(i)*h_grad_eval_new{i};
            end
%             L_grad_eval_new = eval(sym(subs(L,sym_list,sym(xk))));
            
            
            p_quad = xk - xk_prev;
            q_quad = L_grad_quad_new - L_grad_quad;
            
            Ck_quad = (q_quad*q_quad')/(q_quad'*p_quad) - B*p_quad*p_quad'*B/(p_quad'*B*p_quad);
            
            B = B+ Ck_quad;
            
        end
        x=xk; 
                
           

    case 'MSQP'
        % MSQP ALGORITHM HERE
%                           (a) data.x0 = initial solution
%                           (b) data.u0 = inequality constraint Lagrangian
%                           (c) data.v0 = equality constraint Lagrangian
%                           (d) data.B0 = initial Hessian approximation
%                           (e) data.c = convergence parameter
        EPSILON=data.epsilon;
        x=data.x0;
        u=data.u0;
        v=data.v0;
        B=data.B0;
        n=length(x);
        m=length(u);
        l=length(v);
        c=data.c0;
        k=0;
        
        %%% 
        sym_list=[];
        for  i=1:n
            sym_list=[sym_list sym(sprintf('x%d',i))];
        end
        % grad and hessian of f
        f_grad = gradient(data.f(sym_list),sym_list);
        f_hess = hessian(data.f(sym_list),sym_list);

        % grad and hessian of g
        for i=1:m
            g_grad{i} = gradient(data.g{i}(sym_list),sym_list);
            g_hess{i} = hessian(data.g{i}(sym_list),sym_list);
        end
        % grad and hessian of h
        for i=1:l
            h_grad{i} = gradient(data.h{i}(sym_list),sym_list);
            h_hess{i} = hessian(data.h{i}(sym_list),sym_list);            
        end
        

        % to evaluate the gradient and hessian,
        xk=x;
        xk_prev=xk;
        while k<1000
            f_grad_eval=eval(sym(subs(f_grad,sym_list,sym(xk))));
            f_hess_eval=eval(sym(subs(f_hess,sym_list,sym(xk))));
            for i=1:m
                g_grad_eval{i}=eval(sym(subs(g_grad{i},sym_list,sym(xk))));
                g_hess_eval{i}=eval(sym(subs(g_hess{i},sym_list,sym(xk))));
            end
            for i=1:l
                h_grad_eval{i}=eval(sym(subs(h_grad{i},sym_list,sym(xk))));
                h_hess_eval{i}=eval(sym(subs(h_hess{i},sym_list,sym(xk))));
            end
        
            % calculate the H_quad, f_quad
            H_quad=B;
            f_quad=f_grad_eval;
            
            A_quad=[];
            b_quad=[];
            for i=1:m
                A_quad = [A_quad g_grad_eval{i}];
                b_quad = [b_quad; -data.g{i}(xk)];
            end

            
            Aeq_quad=[];
            beq_quad=[];
            for i=1:l
                Aeq_quad = [Aeq_quad h_grad_eval{i}];
                beq_quad = [beq_quad; -data.h{i}(xk)];
            end
            
            %update the hessian matrix
            L_grad_quad = f_grad_eval;
            for i=1:m
                L_grad_quad = L_grad_quad + u(i)*g_grad_eval{i};
            end
            for i=1:l
                L_grad_quad = L_grad_quad + v(i)*h_grad_eval{i};
            end
%             L_grad_eval_prev=eval(sym(subs(L,sym_list,sym(xk))));
            
            [d,fval,exitflag,output,lambda]=quadprog(H_quad,f_quad,A_quad',b_quad,Aeq_quad',beq_quad);
            d;
            if norm(d) < EPSILON
                break;
            end
            xk_prev=xk;
            u = lambda.ineqlin;
            v = lambda.eqlin;            
            
            %%% find the best stepsize using penalty method, % PENALTY ALGORITHM HERE
%             EPSILON2=1e-4;
%             c=1000;
            beta=2;
            q=1;
            %x=data.sym;
            p=@(x) 0;
            for i=1:data.m
                p =@(x) p(x) + max(0,data.g{i}(x))^q;
            end
            for i=1:data.l
                p =@(x) p(x) + abs(data.h{i}(x))^q;
            end
            
            
%             k=0;
%             x=data.x0;
            c = max([data.c0,max(abs(v)),max(u)])
            k = k+1;
%             while c<1e10
            L=@(stepSize) data.f(x + stepSize*d) + c*p(x + stepSize*d); 
            [stepSize,fval]=fminunc(L,1);
%                 if p(x + stepSize*d) < EPSILON
%                     break
%                 end

%             end
            d
%             xk = xk + stepSize * d
            xk = xk +  d;
            fvalAll=[fvalAll data.f(xk)]


            
            % update the hessian matrix using BFGS
            f_grad_eval_new=eval(sym(subs(f_grad,sym_list,sym(xk))));
            for i=1:m
                g_grad_eval_new{i}=eval(sym(subs(g_grad{i},sym_list,sym(xk))));
            end
            for i=1:l
                h_grad_eval_new{i}=eval(sym(subs(h_grad{i},sym_list,sym(xk))));
            end
            
            L_grad_quad_new = f_grad_eval_new;
            for i=1:m
                L_grad_quad_new = L_grad_quad_new + u(i)*g_grad_eval_new{i};
            end
            for i=1:l
                L_grad_quad_new = L_grad_quad_new + v(i)*h_grad_eval_new{i};
            end
%             L_grad_eval_new = eval(sym(subs(L,sym_list,sym(xk))));
            
            
            p_quad = xk - xk_prev;
            q_quad = L_grad_quad_new - L_grad_quad;
            
            Ck_quad = (q_quad*q_quad')/(q_quad'*p_quad) - B*p_quad*p_quad'*B/(p_quad'*B*p_quad);
            
            B = B+ Ck_quad;
            
        end
        x=xk;     
        
    
end

x_optimal=x;
obj_value=data.f(x_optimal);
runtime=toc;
figure
plot(1:length(fvalAll),fvalAll+0.56);
xlabel('Number of iteration') % x-axis label
ylabel('Objective value') % y-axis label
end