mode(1)
//
// Demo of fgoalattain.sci
//

function f1 = gattainObjfun(x)
f1(1)=2*x(1)*x(1)+x(2)*x(2)-48*x(1)-40*x(2)+304
f1(2)=-x(1)*x(1)-3*x(2)*x(2)
f1(3)=x(1)+3*x(2)-18
f1(4)=-x(1)-x(2)
f1(5)=x(1)+x(2)-8
endfunction
x0=[-1,1];
halt()   // Press return to continue
 
goal=[-5,-3,-2,-1,-4];
weight=abs(goal)
gval  =
[- 0.0000011
- 63.999998
- 2.0000002
- 8.
3.485D-08]
z  =
[4.    3.99]
halt()   // Press return to continue
 
Run fgoalattain
[x,fval,attainfactor,exitflag,output,lambda]=fgoalattain(gattainObjfun,x0,goal,weight)
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
