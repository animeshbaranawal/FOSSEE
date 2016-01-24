mode(1)
//
// Demo of fminimaxChecktype.sci
//

// The function takes a string argument.
function myfunction ( x )
fminimaxChecktype ( "myfunction" , x , "x" , 1 , "string" )
disp("This is a string")
endfunction
// Calling sequences which work
myfunction ( "Scilab" )
// Calling sequences which generate an error
myfunction ( 123456 )
halt()   // Press return to continue
 
// The function takes a string or a matrix of doubles argument.
function myfunction ( x )
fminimaxChecktype ( "myfunction" , x , "x" , 1 , [ "string" "constant" ] )
if ( typeof(x) == "string" ) then
disp("This is a matrix of strings")
else
disp("This is a matrix of doubles")
end
endfunction
// Calling sequences which work
myfunction ( "Scilab" )
myfunction ( 123456 )
// Calling sequences which generate an error
myfunction ( uint8(2) )
halt()   // Press return to continue
 
//========= E N D === O F === D E M O =========//
