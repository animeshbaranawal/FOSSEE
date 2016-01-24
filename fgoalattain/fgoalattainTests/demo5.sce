function f5=objfun5(x)
f5(1)=x(1)+x(2)-9*x(1)*x(2);
f5(2)=x(5)*x(1)+x(4)*x(3);
f5(3)=x(5)-6*x(4)*x(3);
f5(4)=-x(1)*x(1)-3*x(2)*x(2)
f5(5)=x(1)+3*x(2)-18
f5(6) = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
f5(7)= x(2)-x(1)*5+x(2)*x(2)
endfunction

x0=[1,3,2,5.9,4];
goal=[9,0,7,7,-5.6,0,3];
weight=[7.7,0,3.2,0,4,6,2];
A=[7.4,0,-3.6,5,2.7];
b=[9];
Aeq=[8,6,-4,0,2];
beq=[3];
lb=[6,3,9,1,2];
ub=[7,4,100,45.9,67.8];

function [c,ceq]=nonlinfun(x)
c=[3*x(4)*x(5)+x(1)];
ceq=[7*x(1)+x(2),x(4)-x(3)*x(3)];
endfunction

[z,gval,attainfactor,exitflag,output,lambda]=fgoalattain(objfun5,x0,goal,weight,A,b,Aeq,beq,lb,ub,nonlinfun)
