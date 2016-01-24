function f1 = objfun(x)
    f1(1)=2*x(1)*x(1)+x(2)*x(2)-48*x(1)-40*x(2)+304
    f1(2)=-x(1)*x(1)-3*x(2)*x(2)
    f1(3)=x(1)+3*x(2)-18
    f1(4)=-x(1)-x(2)
    f1(5)=x(1)+x(2)-8    
endfunction

function G = fgrad(x)
G=[4*x(1)-48,-2*x(1),1,-1,1;2*x(2)-40,-6*x(2),3,-1,1]'
endfunction

x0=[-1,1];
goal=[-5,-3,-2,-1,-4];
weight=abs(goal);
A=[];
b=[];
Aeq=[];
beq=[];
lb=[];
ub=[];

function [c,ceq]=nonlinfun(x)
c=[x(1),-5]
ceq=[x(1)+2*x(2)]
endfunction

function [dc,dceq] = cgrad(x)
dc=[1,0;0,0]';
dceq=[1;2]';
endfunction

options=list("MaxIter", [10000], "CpuTime", [5000], "GradObj", fgrad,"GradCon", cgrad);
[z,gval,attainfactor,exitflag,output,lambda]=fgoalattain(objfun,x0,goal,weight,A,b,Aeq,beq,lb,ub,nonlinfun,options)

