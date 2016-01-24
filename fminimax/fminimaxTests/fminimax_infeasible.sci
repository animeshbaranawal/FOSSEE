function f = myfun(x)
    f(1) = abs(-x(1) - x(2) - x(3))
endfunction

A = [1 2 -1; -1 0 0; 0 -1 0; 0 0 -1]
b = [-4; 0; 0; 0]
Aeq = [1 5 1; 1 1 0]
beq = [10; 100]

x0 = [0 0 0]
[x,fval,maxfval,exitflag,output,lambda] = fminimax(myfun,x0,A,b,Aeq,beq)
