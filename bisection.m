function res = bisection(func, a, b, tol)
    while (b - a) / 2 > tol
        c = (a + b) / 2;
        if func(c) == 0
            break;
        elseif sign(func(c)) == sign(func(a))
            a = c;
        else
            b = c;
        end
    end
%    x_root = (a + b) / 2;
    res = [a,b];
end