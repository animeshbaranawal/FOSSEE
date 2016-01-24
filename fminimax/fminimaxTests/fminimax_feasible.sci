
function f = myfun(x)
    f(1) = 0
    for i = 1:10
        f(1) = f(1) + 4*x(i) 
    end
    for i = 11:20
        f(1) = f(1) - abs(x(i))
    end
endfunction

A = zeros(20,20)
for i = 1:10
    A(i,i) = -1
end
for i = 11:20
    A(i,i) = 1
end

b = zeros(20,1)

function t = tmpf(x)
    t = 0
    for i = 1:10
        t = t+abs(x(i))
    end
    for i = 11:20
        t = t+x(i)
    end
    t = t - 1
endfunction


function [c,ceq] = confun(x)
    c = []
    ceq = tmpf(x)
endfunction

x0 = zeros(1,20)

minimaxoptions = list(..
    "MaxIter"     , [10000], ...
    "CpuTime"   , [10], ...
    "GradObj" ,   "OFF", ...
    "GradCon",   "OFF"...
    );

[x,fval,maxfval,exitflag,output,lambda] = fminimax(myfun,x0,A,b,[],[],[],[],confun,minimaxoptions)
