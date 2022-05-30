function f = obj(x)
    f = (x(1) - 0.2)^2 + (x(2) - 0.5)^2 + (x(3) - 0.3)^2
endfunction

function [varargout, c, ceq] = seminfcon(z,t)
    c = []
    ceq = []
    varargout(1) = sin(t(1)*z(1))*cos(t(1)*z(2)) - 0.001*(t(1)-50)^2 - sin(t(1)*z(3)) - z(3) - 1
    varargout(2) = sin(t(1)*z(2))*cos(t(1)*z(1)) - 0.001*(t(1)-50)^2 - sin(t(1)*z(3)) - z(3) - 1
endfunction

w = [1 0 100 0; 1 0 100 0]
S = [0.2 0; 0.2 0]

[x,fval] = fseminf(obj,[0.5 0.2 0.3],2,seminfcon,w,S,[],[],[],[],[],[])


//function Ki = semiK(z,t,i)
//    [v,t1,t2] = seminfcon(z,t)
//    //disp(i)
//    Ki = v(i)
//endfunction
//
//t = [1:0.2:100]';
//for j = 1:size(t,'r')
//    y(j) = semiK([0.5 0.5 0.5],t(j),1);
//end
//dy = splin(t,y);
//tt = [1:0.2/100:100]';
//yy = interp(tt,t,y,dy);
//Kimax = max(yy);
//
//disp(Kimax);
//
//t = [1:0.2:100]';
//for j = 1:size(t,'r')
//    y(j) = semiK(z,t(j),2);
//end
//dy = splin(t,y);
//tt = [1:0.2/100:100]';
//yy = interp(tt,t,y,dy);
//Kimax = max(yy);
//
//disp(Kimax);
//
