from sympy import symbols, sqrt, exp, ln, pi
from sympy.parsing.sympy_parser import parse_expr


def PW92c_alpha_replaced(parameters):
    # Special version of PW92, used for PKZB functional
    #  PW92(rho_alpha, rho_beta) -> PW92(rho_alpha, 0)
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    t,u,v,w,x,y,p = symbols('t u v w x y p', real=True)
    rho_a = 1/2*(rho_c + rho_s)
    
    # Paramters
    c = parse_expr(parameters["c_PW92c"])#c = 1.709921
    #c = 1.709920934161365617563962776245 # libxc value
    T = [parse_expr(parameters["T1_PW92c"]),parse_expr(parameters["T2_PW92c"]),parse_expr(parameters["T3_PW92c"])]#T = [(1-ln(2))/pi**2, (1-ln(2))/(2*pi**2), 0.016887]# Dalton value
    #T = [0.031091, 0.015545, 0.016887] # paper value
    #T = [0.0310907, 0.01554535, 0.0168869] # libxc value
    U = [parse_expr(parameters["U1_PW92c"]),parse_expr(parameters["U2_PW92c"]),parse_expr(parameters["U3_PW92c"])]#U = [0.21370,  0.20548,  0.11125]
    V = [parse_expr(parameters["V1_PW92c"]),parse_expr(parameters["V2_PW92c"]),parse_expr(parameters["V3_PW92c"])]#V = [7.5957, 14.1189, 10.357]
    W = [parse_expr(parameters["W1_PW92c"]),parse_expr(parameters["W2_PW92c"]),parse_expr(parameters["W3_PW92c"])]#W = [3.5876, 6.1977, 3.6231]
    X = [parse_expr(parameters["X1_PW92c"]),parse_expr(parameters["X2_PW92c"]),parse_expr(parameters["X3_PW92c"])]#X = [1.6382, 3.3662,  0.88026]
    Y = [parse_expr(parameters["Y1_PW92c"]),parse_expr(parameters["Y2_PW92c"]),parse_expr(parameters["Y3_PW92c"])]#Y = [0.49294, 0.62517, 0.49671]
    P = [parse_expr(parameters["P1_PW92c"]),parse_expr(parameters["P2_PW92c"]),parse_expr(parameters["P3_PW92c"])]#P = [1.0, 1.0, 1.0]
    
    # Functional
    zeta = rho_a/rho_a
    r = (3/(4*pi*rho_a))**(1/3)
    omega =((1+zeta)**(4/3)+(1-zeta)**(4/3)-2)/(2**(4/3)-2)
    en = -2*t*(1+u*r)*ln(1+1/(2*t*(v*sqrt(r)+w*r+x*r**(3/2)+y*r**(p+1))))
    eps = en.subs({t:T[0],u:U[0],v:V[0],w:W[0],x:X[0],y:Y[0],p:P[0]})
    eps = eps - en.subs({t:T[2],u:U[2],v:V[2],w:W[2],x:X[2],y:Y[2],p:P[2]})*omega*(1-zeta**4)/c
    eps = eps + (en.subs({t:T[1],u:U[1],v:V[1],w:W[1],x:X[1],y:Y[1],p:P[1]}) -en.subs({t:T[0],u:U[0],v:V[0],w:W[0],x:X[0],y:Y[0],p:P[0]}))*omega*zeta**4
    
    return eps
    
    
def PW92c_beta_replaced(parameters):
    # Special version of PW92, used for PKZB functional
    #  PW92(rho_alpha, rho_beta) -> PW92(rho_beta, 0)
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    t,u,v,w,x,y,p = symbols('t u v w x y p', real=True)
    rho_b = 1/2*(rho_c - rho_s)
    
    # Paramters
    c = parse_expr(parameters["c_PW92c"])#c = 1.709921
    #c = 1.709920934161365617563962776245 # libxc value
    T = [parse_expr(parameters["T1_PW92c"]),parse_expr(parameters["T2_PW92c"]),parse_expr(parameters["T3_PW92c"])]#T = [(1-ln(2))/pi**2, (1-ln(2))/(2*pi**2), 0.016887]# Dalton value
    #T = [0.031091, 0.015545, 0.016887] # paper value
    #T = [0.0310907, 0.01554535, 0.0168869] # libxc value
    U = [parse_expr(parameters["U1_PW92c"]),parse_expr(parameters["U2_PW92c"]),parse_expr(parameters["U3_PW92c"])]#U = [0.21370,  0.20548,  0.11125]
    V = [parse_expr(parameters["V1_PW92c"]),parse_expr(parameters["V2_PW92c"]),parse_expr(parameters["V3_PW92c"])]#V = [7.5957, 14.1189, 10.357]
    W = [parse_expr(parameters["W1_PW92c"]),parse_expr(parameters["W2_PW92c"]),parse_expr(parameters["W3_PW92c"])]#W = [3.5876, 6.1977, 3.6231]
    X = [parse_expr(parameters["X1_PW92c"]),parse_expr(parameters["X2_PW92c"]),parse_expr(parameters["X3_PW92c"])]#X = [1.6382, 3.3662,  0.88026]
    Y = [parse_expr(parameters["Y1_PW92c"]),parse_expr(parameters["Y2_PW92c"]),parse_expr(parameters["Y3_PW92c"])]#Y = [0.49294, 0.62517, 0.49671]
    P = [parse_expr(parameters["P1_PW92c"]),parse_expr(parameters["P2_PW92c"]),parse_expr(parameters["P3_PW92c"])]#P = [1.0, 1.0, 1.0]
    
    # Functional
    zeta = rho_b/rho_b
    r = (3/(4*pi*rho_b))**(1/3)
    omega =((1+zeta)**(4/3)+(1-zeta)**(4/3)-2)/(2**(4/3)-2)
    en = -2*t*(1+u*r)*ln(1+1/(2*t*(v*sqrt(r)+w*r+x*r**(3/2)+y*r**(p+1))))
    eps = en.subs({t:T[0],u:U[0],v:V[0],w:W[0],x:X[0],y:Y[0],p:P[0]})
    eps = eps - en.subs({t:T[2],u:U[2],v:V[2],w:W[2],x:X[2],y:Y[2],p:P[2]})*omega*(1-zeta**4)/c
    eps = eps + (en.subs({t:T[1],u:U[1],v:V[1],w:W[1],x:X[1],y:Y[1],p:P[1]}) -en.subs({t:T[0],u:U[0],v:V[0],w:W[0],x:X[0],y:Y[0],p:P[0]}))*omega*zeta**4
    
    return eps
    
    
def PBEc_alpha_replaced(parameters):
    # Special version of PBE, used for PKZB functional
    #  PBE(rho_alpha, rho_beta, grad_alpha, grad_beta) -> PBE(rho_alpha, 0, grad_alpha, 0)
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    gamma_cc, gamma_ss = symbols('gamma_cc gamma_ss', real=True, nonnegative=True)
    gamma_cs = symbols('gamma_cs', real=True)
    rho_a = 1/2*(rho_c + rho_s)
    gamma_aa = 1/4*(gamma_cc + gamma_ss + 2*gamma_cs)
    
    # Parameters
    gamma = parse_expr(parameters["gamma_PBEc"])#gamma = (1 - ln(2))/(pi**2)
    beta = parse_expr(parameters["beta_PBEc"])#beta = 0.066725
    #beta = 0.06672455060314922 # libxc value
    
    # Functional
    a0 = 1
    e = 1
    kF = (3 * rho_a * pi**2)**(1/3)
    ks = sqrt( 4 * kF/(pi*a0) )
    zeta = rho_a/rho_a
    phi = ( (1 + zeta)**(2/3) + (1 - zeta)**(2/3) )/2
    t = sqrt(gamma_aa)/( 2 * phi * ks * rho_a )
    t2 = t**2
    ecunif = PW92c_alpha_replaced(parameters)
    A = beta/gamma*(exp(-ecunif/(gamma*phi**3*e**2/a0))-1)**(-1)
    H = (e**2/a0)*gamma*phi**3*ln(1 + beta/gamma*t**2*(( 1 + A * t**2 )/( 1 + A * t**2 + A**2 * t**4 )))
    epsPBE = ecunif + H
    
    return epsPBE
    
    
def PBEc_beta_replaced(parameters):
    # Special version of PBE, used for PKZB functional
    #  PBE(rho_alpha, rho_beta, grad_alpha, grad_beta) -> PBE(rho_beta, 0, grad_beta, 0)
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    gamma_cc, gamma_ss = symbols('gamma_cc gamma_ss', real=True, nonnegative=True)
    gamma_cs = symbols('gamma_cs', real=True)
    rho_b = 1/2*(rho_c - rho_s)
    gamma_bb = 1/4*(gamma_cc + gamma_ss - 2*gamma_cs)
    
    # Parameters
    gamma = parse_expr(parameters["gamma_PBEc"])#gamma = (1 - ln(2))/(pi**2)
    beta = parse_expr(parameters["beta_PBEc"])#beta = 0.066725
    #beta = 0.06672455060314922 # libxc value
    
    # Functional
    a0 = 1
    e = 1
    kF = (3 * rho_b * pi**2)**(1/3)
    ks = sqrt( 4 * kF/(pi*a0) )
    zeta = rho_b/rho_b
    phi = ( (1 + zeta)**(2/3) + (1 - zeta)**(2/3) )/2
    t = sqrt(gamma_bb)/( 2 * phi * ks * rho_b )
    t2 = t**2
    ecunif = PW92c_beta_replaced(parameters)
    A = beta/gamma*(exp(-ecunif/(gamma*phi**3*e**2/a0))-1)**(-1)
    H = (e**2/a0)*gamma*phi**3*ln(1 + beta/gamma*t**2*(( 1 + A * t**2 )/( 1 + A * t**2 + A**2 * t**4 )))
    epsPBE = ecunif + H
    
    return epsPBE
