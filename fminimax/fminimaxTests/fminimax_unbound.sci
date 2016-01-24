
function f = myfun(x)
    f(1) = x(1)^3 - 2
    f(2) = x(1)^3 
    f(3) = x(1)^3 + 6
endfunction

x0 = [2]
[x,fval,maxfval,exitflag,output,lambda] = fminimax(myfun,x0)