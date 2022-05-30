function [x,fval] = fseminf(fun,x0,ntheta,seminfcon,w,s,A,b,Aeq,beq,lb,ub)

    fsinfObjfun = fun
    fsinfStartpoint = x0
    intervals = w
    samplespace = s
    fsinfA = A
    fsinfb = b
    fsinfAeq = Aeq
    fsinfBeq = beq
    fsinflb = lb
    fsinfub = ub

    [v,t1,t2] = seminfcon(fsinfStartpoint,[0 0])
    //disp(v)
    //disp(size(v,'c')) 
    //disp(size(v,'r'))
    assert_checktrue(size(v,'r') == ntheta)
    assert_checktrue(size(w,'r') == ntheta)
    assert_checktrue(size(w,'c') == 4)
    assert_checktrue(size(s,'r') == ntheta)
    assert_checktrue(size(s,'c') == 2)

    function Ki = semiK(z,t,i)
        [v,t1,t2] = seminfcon(z,t)
        Ki = v(i)
    endfunction

    function Kimax = findMaxi(z,i)
        t = [intervals(i,1):samplespace(i,1):intervals(i,3)]';
        for j = 1:size(t,'r')
            y(j) = semiK(z,t(j),i);
        end
        dy = splin(t,y);
        tt = [intervals(i,1):samplespace(i,1)/100:intervals(i,3)]';
        yy = interp(tt,t,y,dy);
        Kimax = max(yy);
    endfunction

    function f = allMax(z)
        f = [];
        for i = 1:ntheta
            //disp(findMaxi(z,i));
            f = [f; findMaxi(z,i)];
        end
        f = max(f);
    endfunction
    
    funcprot(0)
    isMaxPos = 0
    isMaxNeg = 0
    MaxIter = 100

    xStart = fsinfStartpoint
    interpmax = allMax(xStart)
    
    disp(xStart)
    disp(interpmax)
    
    function [c,ceq] = fsinfNonlinfun(z)
        [t1,c,ceq] = seminfcon(z,[0 0]);
        tmpVal = allMax(z)
        if interpmax > 0 then
            c = [c, interpmax];
        end
    endfunction
    xNext = fmincon(fsinfObjfun,xStart,fsinfA,fsinfb,fsinfAeq,fsinfBeq,fsinflb,fsinfub,fsinfNonlinfun,list("MaxIter",1))

    if(allMax(xNext) > 0) then
        isMaxPos = 1
    else
        isMaxNeg = 1
    end

    exitFlag = 0

    if isMaxNeg == 1 then

        while(exitFlag == 0 & MaxIter > 0)
            disp([1 2 3 4])
            xStart = xNext
            interpmax = allMax(xStart)
            
            disp(xStart)
            disp(interpmax)
            
            function [c,ceq] = fsinfNonlinfun(z)
                [t1,c,ceq] = seminfcon(z,[0 0]);
                tmpVal = allMax(z);
                if interpmax > 0 then
                    c = [c, interpmax];
                end
            endfunction
            xNext = fmincon(fsinfObjfun,xStart,fsinfA,fsinfb,fsinfAeq,fsinfBeq,fsinflb,fsinfub,fsinfNonlinfun,list("MaxIter",1))

            MaxIter = MaxIter - 1
            if allMax(xNext) > 0 then
                disp([1 2 3 4 5 6 7 8 9])
                exitFlag = 1
            end
        end

    else

        while(exitFlag == 0 & MaxIter > 0)
            disp([1 2 3 5])
            xStart = xNext
            interpmax = allMax(xStart)
            
            disp(xStart)
            disp(interpmax)
            
            function [c,ceq] = fsinfNonlinfun(z)
                [t1,c,ceq] = seminfcon(z,[0 0]);
                tmpVal = allMax(z)
                if interpmax > 0 then
                    c = [c, interpmax];
                end
            endfunction
            xNext = fmincon(fsinfObjfun,xStart,fsinfA,fsinfb,fsinfAeq,fsinfBeq,fsinflb,fsinfub,fsinfNonlinfun,list("MaxIter",1))

            MaxIter = MaxIter - 1
            if allMax(xNext) <= 0 then
                disp([1 2 3 4 5 6 7 8 9])
                exitFlag = 1
            end
        end

        exitFlag = 0
        MaxIter = 10000
        while(exitFlag == 0 & MaxIter > 0)
            disp([1 2 3 6])
            xStart = xNext
            interpmax = allMax(xStart)
            
            disp(xStart)
            disp(interpmax)
            
            function [c,ceq] = fsinfNonlinfun(z)
                [t1,c,ceq] = seminfcon(z,[0 0]);
                tmpVal = allMax(z)
                if interpmax > 0 then
                    c = [c, interpmax];
                end
            endfunction
            xNext = fmincon(fsinfObjfun,xStart,fsinfA,fsinfb,fsinfAeq,fsinfBeq,fsinflb,fsinfub,fsinfNonlinfun,list("MaxIter",1))

            MaxIter = MaxIter - 1
            if allMax(xNext) > 0 then
                exitFlag = 1
                disp([1 2 3 4 5 6 7 8 9])
            end
        end

    end

    fval = fsinfObjfun(xStart)
    x = xStart

endfunction
