function f = obj(x)
    f = (x(1) - 2)^2
endfunction

function [varargout,c,ceq] = seminfcon(x,t)
    c = []
    ceq = []
    varargout(1) = (x(1) - 0.5) - (t(1) - 0.5)^2
endfunction

[x,fval] = fseminf(obj,[3],1,seminfcon,[0 0 1 0],[0.0001 0],[],[],[],[],[],[])
