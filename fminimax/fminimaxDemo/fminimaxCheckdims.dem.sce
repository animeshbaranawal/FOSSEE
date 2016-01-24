mode(1)
//
// Demo of fminimaxCheckdims.sci
//

// The function takes a 2 x 3 matrix of doubles.
function y = myfunction ( x )
fminimaxCheckdims ( "myfunction" , x , "x" , 1 , [2 3] )
y = x
endfunction
// Calling sequences which work
y = myfunction ( ones(2,3) )
y = myfunction ( zeros(2,3) )
// Calling sequences which generate an error
y = myfunction ( ones(1,3) )
y = myfunction ( zeros(2,4) )
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
