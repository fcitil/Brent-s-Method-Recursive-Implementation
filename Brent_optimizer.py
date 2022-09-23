from ast import main
import numpy as np


lookup = dict() 
max_iter = 100
eps = 1e-6
InvPhi =(1 + np.sqrt(5)) / 2 - 1
ratio = (3 - np.sqrt(5)) / 2


def f(x,fun):
    """
    This function is used to store and use the values of the function in the form of (x,fun): fun(x)
    param x: value of x
    param fun: function to minimize
    return: the value of the function
    """
    if (x,fun) in lookup:
        return lookup[(x,fun)]
    else:
        lookup[(x,fun)] = fun(x)
        return lookup[(x,fun)]

# def jarratt(xnm2, xnm1, xn, fun):
#     """
#     Tis funciton uses successive parabolic interpolation to find the minimum of the function
#     param xnm2: x value at n-2
#     param xnm1: x value at n-1
#     param xn: x value at n
#     param fun: function to minimize
#     return: the minimum of the function
#     """
#     fn = f(xn,fun)
#     fnm1 = f(xnm1,fun)
#     fnm2 = f(xnm2,fun)
#     return xn + 0.5 * ((xnm1 - xn)**2.0 * (fn - fnm2) + (xnm2 - xn)**2.0 * (fnm1 - fn)) / ((xnm1 - xn) * (fn - fnm2) + (xnm2 - xn) * (fnm1 - fn))

# def jarrattList(x,fun):
#     """
#     One iteration of the Jarratt method using a list
#     param x: list of x values
#     param fun: function to minimize
#     return: the result of a single iteration of the Jarratt method
#     """
#     n = len(x)
#     x.append(jarratt(x[n-3], x[n-2], x[n-1], fun))
#     return x
    
# def polynomialInterpolation(x, fun):
#     """
#     Jarratt's Method iterarting until xmin is found
#     param x: list of x values
#     param fun: function to minimize
#     return: the minimum of the function
#     """
#     n = len(x)
#     if abs(x[n-1]-x[n-2]) < eps or n>max_iter:
#         return x
#     else:
#         return polynomialInterpolation(jarrattList(x,fun),fun)

# def golden(a,c,d,b,fun):
#     """
#     Golden section search method
#     param a: left bound
#     param c: first test point
#     param d: second test point
#     param b: right bound
#     return: the minimum of the function
#     """
#     if b - a < eps:
#         return (a + b) / 2
    
#     elif fun(c) < fun(d):
#         t = d + InvPhi*(a-d)
#         return golden(a,t,c,d,fun)

#     else:
#         t = c + InvPhi*(b-c)
#         return golden(c,d,t,b,fun)
# def golden_Section(a, b, fun):
#     """
#     Golden-section Search Initializer
#     param a: left bound
#     param b: right bound
#     return: the minimum of the function
#     """
#     c = b + InvPhi * (a - b)
#     d = a + InvPhi * (b - a)

#     return golden(a,c,d,b,fun)


def brent(a, b, v, w, x, d_old, e_old, i, fun):
    """
    Brent's Method
    param a: left bound
    param b: right bound
    param v: previous iterate value
    param w: previous iterate value
    param x: previous iterate value
    param d_old: last delta step
    param e_old: last golden interval size
    param i: iteration counter
    param fun: function to minimize
    return: the minimum of the function
    """
    fv = f(v,fun)
    fw = f(w,fun)
    fx = f(x,fun)
    new_i = i + 1
    m = 0.5 * (a + b)
    print("Iteration: ", new_i)
    print("epsilon: ",(b-a)/2)
    if b-a <= eps or i > max_iter:
        print("Minimum found at: ", x)
        return m
    
    else:
        r = (x - w) * (fx - fv)
        tq = (x - v) * (fx - fw)
        tp = (x - v) * tq - (x - w) * r
        tq_2 = 2 * (tq - r)
        p = -tp if tq_2 > 0 else tp
        q = tq_2 if tq_2 > 0 else -tq_2
        safe = q != 0
        delta_x = p/q if safe else 0
        parabolic = safe and a < (x + delta_x) and x + delta_x < b and abs(delta_x) < 0.5 * abs(e_old)
        
        if parabolic:
            print("Parabolic")
            e = d_old
        elif x < m:
            print("Golden")
            e = b - x
        else:
            print("Golden_1")
            e = a - x
        
        d = delta_x if parabolic else ratio * e
        u = x + d
        fu = f(u,fun)
        if fu <= fx:
            if u < x:
                new_a = a
                new_b = x
            else:
                new_a = x
                new_b = b

            return brent(new_a, new_b, w, x, u, d, e, new_i, fun)
        else:
            if u < x:
                new_a = u
                new_b = x

            else:
                new_a = a
                new_b = u  

            if fu <= fx or w == x:
                return brent(new_a, new_b, w, u, x, d, e, new_i, fun)
            
            elif fu <= fw or v == x or v == w:
                return brent(new_a, new_b, u, w, x, d, e, new_i, fun)
            
            else: 
                return brent(new_a, new_b, v, w, x, d, e, new_i, fun)

def brent_Method(fun, bounds, xtol=1e-6, maxiter=100):
    """
    Brent's Method Initializer
    param a: left bound
    param b: right bound
    param fun: function to minimize
    return: the minimum of the function
    """
    global eps, max_iter

    eps = xtol
    max_iter = maxiter

    a, b = bounds
    x = b + InvPhi * (a - b)
    return brent(a, b, x, x, x, 0, 0, 0, fun)


def test_function(x):
    """
    Test Function to be minimized
    """
    return  x**4 / 4 - x**3 / 2 - x**2 / 2 - 1 



if __name__ == "__main__":
    print("Brent's Method: ", brent_Method(test_function, (1, 2)))