mode(1)
//
// Demo of fminimaxCheckrhs.sci
//

// The function takes 2/3 input arguments and 1 output arguments
function y = myfunction ( varargin )
[lhs, rhs] = argn()
fminimaxCheckrhs ( "myfunction" , rhs , 2:3 )
fminimaxChecklhs ( "myfunction" , lhs , 1 )
x1 = varargin(1)
x2 = varargin(2)
if ( rhs >= 3 ) then
x3 = varargin(3)
else
x3 = 2
end
y = x1 + x2 + x3
endfunction
// Calling sequences which work
y = myfunction ( 1 , 2 )
y = myfunction ( 1 , 2 , 3 )
// Calling sequences which generate an error
y = myfunction ( 1 )
y = myfunction ( 1 , 2 , 3 , 4 )
halt()   // Press return to continue
 
// The function takes 2 or 4 input arguments, but not 3
function y = myfunction ( varargin )
[lhs, rhs] = argn()
fminimaxCheckrhs ( "myfunction" , rhs , [2 4] )
fminimaxChecklhs ( "myfunction" , lhs , 1 )
x1 = varargin(1)
x2 = varargin(2)
if ( rhs >= 3 ) then
x3 = varargin(3)
x4 = varargin(4)
else
x3 = 2
x4 = 3
end
y = x1 + x2 + x3 + x4
endfunction
// Calling sequences which work
y = myfunction ( 1 , 2 )
y = myfunction ( 1 , 2 , 3 , 4 )
// Calling sequences which generate an error
y = myfunction ( 1 )
y = myfunction ( 1 , 2 , 3 )
y = myfunction ( 1 , 2 , 3 , 4, 5 )
halt()   // Press return to continue
 
// The function takes 2 input arguments and 0/1 output arguments.
// Notice that if the checkrhs function is not called,
// the variable x2 might be used from the user's context,
// that is, if the caller has defined the variable x2, it
// is used in the callee.
// Here, we want to avoid this.
function y = myfunction ( x1 , x2 )
[lhs, rhs] = argn()
fminimaxCheckrhs ( "myfunction" , rhs , 2 )
fminimaxChecklhs ( "myfunction" , lhs , [0 1] )
y = x1 + x2
endfunction
// Calling sequences which work
y = myfunction ( 1 , 2 )
// Calling sequences which generate an error
y = myfunction ( 1 )
y = myfunction ( 1 , 2 , 3 )
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
