mode(1)
//
// Demo of fminimaxChecklhs.sci
//

// The function takes 3 input arguments and 1/2 output arguments
function varargout = myfunction ( x1 , x2 , x3 )
[lhs, rhs] = argn()
fminimaxCheckrhs ( "myfunction" , rhs , 3 : 3 )
fminimaxChecklhs ( "myfunction" , lhs , 1 : 2 )
y1 = x1 + x2
y2 = x2 + x3
varargout(1) = y1
if ( lhs == 2 ) then
varargout(2) = y2
end
endfunction
// Calling sequences which work
myfunction ( 1 , 2 , 3 )
y1 = myfunction ( 1 , 2 , 3 )
[ y1 , y2 ] = myfunction ( 1 , 2 , 3 )
// Calling sequences which generate an error
[ y1 , y2 , y3 ] = myfunction ( 1 , 2 , 3 )
halt()   // Press return to continue
 
// The function takes 1 or 3 output arguments, but not 2
function varargout = myfunction ( x1 , x2 , x3 )
[lhs, rhs] = argn()
fminimaxCheckrhs ( "myfunction" , rhs , 3 : 3 )
fminimaxChecklhs ( "myfunction" , lhs , [1 3] )
y1 = x1 + x2
y2 = x2 + x3
y3 = x1 + x3
varargout(1) = y1
if ( lhs == 3 ) then
varargout(2) = y2
varargout(3) = y3
end
endfunction
// Calling sequences which work
myfunction ( 1 , 2 , 3 )
y1 = myfunction ( 1 , 2 , 3 )
[ y1 , y2 , y3 ] = myfunction ( 1 , 2 , 3 )
// Calling sequences which generate an error
[y1 , y2] = myfunction ( 1 , 2 , 3 )
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
