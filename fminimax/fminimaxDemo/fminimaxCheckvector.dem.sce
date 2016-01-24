mode(1)
//
// Demo of fminimaxCheckvector.sci
//

// The function takes a vector of 3 doubles.
function y = myfunction ( x )
fminimaxCheckvector ( "myfunction" , x , "x" , 1 , 3 )
y = x
endfunction
// Calling sequences which work
y = myfunction ( ones(1,3) )
y = myfunction ( zeros(3,1) )
// Calling sequences which generate an error
// The following are not vectors
y = myfunction ( ones(2,3) )
y = myfunction ( zeros(3,2) )
// The following have the wrong number of entries
y = myfunction ( ones(1,3) )
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
