%MySolve
function [x,converged,J] = MySolve(f,x0,df,tol,maxit) 
    x=x0;
    converged=0; %Default value of convergence flag
    for i = 1:maxit
        x_prev = x;
        J=df(x);
        x = x_prev - J\f(x);
        %Runs for maxit iterations or until relative and
        %absolute convergence conditions are met
        if norm(f(x), Inf) <= tol && norm(x-x_prev, Inf) <= tol 
            disp('converged')
            converged = 1;
            break
        else
            disp('iteration');
        end
    end