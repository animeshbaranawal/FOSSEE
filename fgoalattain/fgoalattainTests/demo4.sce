function f4=objfun4(x)
f4(1)=x(2)-x(4)*x(4)*x(3);
f4(2)=x(1)*x(3)-x(2)*x(2);
f4(3)=x(5)*x(1)+x(4)*x(3);
f4(4)=x(3)^2+x(1)-x(2)*x(3);
endfunction

x0=[3,1,-8,-3,0];
goal=[9,0,7,9];
weight=[5,8,0,8];
A=[7,0,-3,5,2];
b=[6];
Aeq=[8,6,-4,0,2];
beq=[9];
lb=[4,6,1,7,6];
ub=[10,11,12,13,14];

[z,gval,attainfactor,exitflag,output,lambda]=fgoalattain(objfun4,x0,goal,weight,A,b,Aeq,beq,lb,ub)