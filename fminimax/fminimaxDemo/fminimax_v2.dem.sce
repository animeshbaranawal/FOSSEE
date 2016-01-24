mode(1)
//
// Demo of fminimax_v2.sci
//

// A basic case :
// we provide only the objective function and the nonlinear constraint
// function
function f = myfun(x)
f(1)= 2*x(1)^2 + x(2)^2 - 48*x(1) - 40*x(2) + 304;     //Objectives
f(2)= -x(1)^2 - 3*x(2)^2;
f(3)= x(1) + 3*x(2) -18;
f(4)= -x(1) - x(2);
f(5)= x(1) + x(2) - 8;
endfunction
halt()   // Press return to continue
 
// The initial guess
x0 = [0.1,0.1];
// The expected solution : only 4 digits are guaranteed
xopt = [4 4]
fopt = [0 -64 -2 -8 0]
maxfopt = 0
// Run fminimax
[x,fval,maxfval,exitflag,output,lambda] = fminimax(myfun, x0)
halt()   // Press return to continue
 
// A case where we provide the gradient of the objective
// functions and the Jacobian matrix of the constraints.
// The objective function and its gradient
function f = myfun(x)
f(1)= 2*x(1)^2 + x(2)^2 - 48*x(1) - 40*x(2) + 304;
f(2)= -x(1)^2 - 3*x(2)^2;
f(3)= x(1) + 3*x(2) -18;
f(4)= -x(1) - x(2);
f(5)= x(1) + x(2) - 8;
endfunction
halt()   // Press return to continue
 
// Defining gradient of myfun
function G = myfungrad(x)
G = [ 4*x(1) - 48, -2*x(1), 1, -1, 1;
2*x(2) - 40, -6*x(2), 3, -1, 1; ]'
endfunction
halt()   // Press return to continue
 
halt()   // Press return to continue
 
// The nonlinear constraints and the Jacobian
// matrix of the constraints
function [c,ceq] = confungrad(x)
// Inequality constraints
c(1) = 1.5 + x(1)*x(2) - x(1) - x(2)
c(2) = -x(1)*x(2) - 10
// No nonlinear equality constraints
ceq=[]
endfunction
halt()   // Press return to continue
 
halt()   // Press return to continue
 
// Defining gradient of confungrad
function [DC,DCeq] = cgrad(x)
// DC(:,i) = gradient of the i-th constraint
// DC = [
//   Dc1/Dx1  Dc1/Dx2
//   Dc2/Dx1  Dc2/Dx2
//   ]
DC= [
x(2)-1, -x(2)
x(1)-1, -x(1)
]'
DCeq = []'
endfunction
halt()   // Press return to continue
 
// Test with both gradient of objective and gradient of constraints
minimaxOptions = list("GradObj","ON","GradCon","ON");
// The initial guess
x0 = [0,10];
// The expected solution : only 4 digits are guaranteed
xopt = [0.92791 7.93551]
fopt = [6.73443  -189.778  6.73443  -8.86342  0.86342]
maxfopt = 6.73443
// Run fminimax
[x,fval,maxfval,exitflag,output] = fminimax(myfun,x0,[],[],[],[],[],[], confungrad,minimaxOptions,myfungrad,cgrad)
halt()   // Press return to continue
 
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
