clear
solver_method='alternate_broyden'; % - 'steepest_descent' - 'newtons' - 'broyden' - 'alternate_broyden'

% %%% small problem
% data.objective=@objective_st;
% data.deri1=@deri1_st;
% data.deri2=@deri2_st;
% data.x0=[0;0]

% %%%large problem
data.objective=@objective_lg;
data.deri1=@deri1_lg;
data.deri2=@deri2_lg;
data.x0=200*(rand(100,1)-0.5); 


data.e=1e-6;
data.e1=1e-3;
data.e2=1e-3;
data.a=1e-10;
data.b=10;
data.phi=0;
data.B0=eye(length(data.x0));
data.D0=eye(length(data.x0));
data.alpha=2;
data.lambda0=0.2;
data.epsilon=0.2;
data.line_search='dich'; % possible choice is 'armijo', 'dich', 'golden','fibonacci','bisection','newton_search'

%%%note; if using newton's method, the algorithm may not converge if the
%%%intitial point is far away from the optimal
% x_optimal = uo_solver2(solver_method, data);
% method=['steepest_descent ';'broyden          ';'alternate_broyden']; %length=17
% linesearch=['armijo       '; 'dich         '; 'golden       ';'fibonacci    ';'bisection    ';'newton_search'];% length=13
% c1=cellstr(method);
% c2=cellstr(linesearch)
% 
% for i=1:3
%     solver_method=c1{i};
%     for j=1:6
%         data.line_search=c2{j};
%         x_optimal = uo_solver(solver_method, data);
%         disp(solver_method),disp('+'),disp(data.line_search),disp('finished')
%         pause(1)
%         
%     end
% end
solver_method='alternate_broyden';
data.line_search='dich'
uo_solver(solver_method, data); %%% mind the initial point x0, if it is two far away from the optimal, then it will not converge