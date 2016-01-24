// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// Authors: Prajwala TM,Sheetal Shalini
// Organization: FOSSEE, IIT Bombay
// Email: prajwala.tm@gmail.com,sheetalsh456@gmail.com
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

// Define subsidiary functions used by fgoalattain

// to check if the number of input arguments is equal to the numbers mentioned in the set
function errmsg = fgoalattainCheckrhs ( funname , rhs , rhsset )  
 errmsg = []
  if ( and(rhs <> rhsset) ) then
    rhsstr = strcat(string(rhsset)," ")
    errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while the number of expected input arguments should be in the set [%s]."), funname , rhs , rhsstr );
    error(errmsg)
  end
endfunction

// to check if the number of output arguments is equal to the numbers mentioned in the set, which is any number between 1 to 6 (both inclusive)
function errmsg = fgoalattainChecklhs ( funname , lhs , lhsset )
errmsg = []
  if ( and ( lhs <> lhsset ) ) then
    lhsstr = strcat(string(lhsset)," ")
    errmsg = msprintf(gettext("%s: Unexpected number of output arguments : %d provided while the expected number of output arguments should be in the set [%s]."), funname , lhs , lhsstr );
    error(errmsg)
  end
endfunction

// Generates an error if the given variable is not of expected type.
function errmsg = fgoalattainChecktype ( funname , var , varname , ivar , expectedtype )
errmsg = []
  if ( and ( typeof ( var ) <> expectedtype ) ) then
    strexp = """" + strcat(expectedtype,""" or """) + """"
    errmsg = msprintf(gettext("%s: Expected type [%s] for input argument %s at input #%d, but got ""%s"" instead."),funname, strexp, varname , ivar , typeof(var) );
    error(errmsg);
  end
endfunction

// Generates an error if the variable is not a vector.
function errmsg = fgoalattainCheckvector ( funname , var , varname , ivar , nbval )
errmsg = []
  nrows = size(var,"r")
  ncols = size(var,"c")
  if ( nrows <> 1 & ncols <> 1 ) then
    strcomp = strcat(string(size(var))," ")
    errmsg = msprintf(gettext("%s: Expected a vector matrix for input argument %s at input #%d, but got [%s] instead."), funname, varname , ivar , strcomp );
    error(errmsg)
  end
  if ( ( nrows == 1 & ncols <> nbval ) | ( ncols == 1 & nrows <> nbval ) ) then
    strcomp = strcat(string(size(var))," ")
    errmsg = msprintf(gettext("%s: Expected %d entries for input argument %s at input #%d, but current dimensions are [%s] instead."), funname, nbval , varname , ivar , strcomp );
    error(errmsg)
  end
endfunction

// Generates an error if the variable has not the required size.
function errmsg = fgoalattainCheckdims ( funname , var , varname , ivar , matdims )
[lhs,rhs]=argn()
  fgoalattainCheckrhs ( "fgoalattainCheckdims" , rhs , 5 )
  fgoalattainChecklhs ( "fgoalattainCheckdims" , lhs , [0 1] )
  errmsg = []
  if ( or ( size(var) <> matdims ) ) then
    strexp = strcat(string(matdims)," ")
    strcomp = strcat(string(size(var))," ")
    errmsg = msprintf(gettext("%s: Expected size [%s] for input argument %s at input #%d, but got [%s] instead."), funname, strexp, varname , ivar , strcomp );
    error(errmsg)
  end
endfunction


