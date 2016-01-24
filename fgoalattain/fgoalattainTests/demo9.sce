function f1 = objfun(x)
    f1(1)=2*x(1)*x(1)+x(2)*x(2)-48*x(1)-40*x(2)+304
    f1(2)=-x(1)*x(1)-3*x(2)*x(2)
    f1(3)=x(1)+3*x(2)-18
    f1(4)=-x(1)-x(2)-1
    f1(5)=x(1)+x(2)-8
    f1(6)=100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
    f1(7)=x(2)-x(1)*5+x(2)*x(2)
    f1(8)=x(2)-4*x(1)+3
    f1(9)=x(2)^3-x(1)^2+4
endfunction

goal=[5.0,-6.3,7.4,2.1,4,6,5.2,8,-1.4]
weight=[8,2,3,4,4,5,3,7,6]
x0=[-1,2]
A=[1,2]
b=[3]
Aeq=[-1,4]
beq=[5]
lb=[-1,-1]
ub=[10,10]

function [p,q]=nonlinfun(x)
p=[3*x(2)^2-2*x(1)];
q=[4*x(2)-7*x(1)];
endfunction

[z,gval,attainfactor,exitflag,output,lambda]=fgoalattain(objfun,x0,goal,weight,A,b,Aeq,beq,lb,ub,nonlinfun)
