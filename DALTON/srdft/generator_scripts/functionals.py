from sympy import symbols, sqrt, exp, ln, pi, atan
from sympy.parsing.sympy_parser import parse_expr
import functionals_special as spec_func


def PW92c(parameters):
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    t,u,v,w,x,y,p = symbols('t u v w x y p', real=True)
    
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
    zeta = rho_s/rho_c
    r = (3/(4*pi*rho_c))**(1/3)
    omega =((1+zeta)**(4/3)+(1-zeta)**(4/3)-2)/(2**(4/3)-2)
    en = -2*t*(1+u*r)*ln(1+1/(2*t*(v*sqrt(r)+w*r+x*r**(3/2)+y*r**(p+1))))
    eps = en.subs({t:T[0],u:U[0],v:V[0],w:W[0],x:X[0],y:Y[0],p:P[0]})
    eps = eps - en.subs({t:T[2],u:U[2],v:V[2],w:W[2],x:X[2],y:Y[2],p:P[2]})*omega*(1-zeta**4)/c
    eps = eps + (en.subs({t:T[1],u:U[1],v:V[1],w:W[1],x:X[1],y:Y[1],p:P[1]}) -en.subs({t:T[0],u:U[0],v:V[0],w:W[0],x:X[0],y:Y[0],p:P[0]}))*omega*zeta**4
    
    return eps
    
    
def VWN5c(parameters):
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    A, b, c, x0 = symbols('A b c x0', real=True)
    
    A1_vwn = parse_expr(parameters["A1_vwn"])
    A2_vwn = parse_expr(parameters["A2_vwn"])
    A3_vwn = parse_expr(parameters["A3_vwn"])
    B1_vwn = parse_expr(parameters["B1_vwn"])
    B2_vwn = parse_expr(parameters["B2_vwn"])
    B3_vwn = parse_expr(parameters["B3_vwn"])
    c1_vwn = parse_expr(parameters["c1_vwn"])
    c2_vwn = parse_expr(parameters["c2_vwn"])
    c3_vwn = parse_expr(parameters["c3_vwn"])
    x01_vwn = parse_expr(parameters["x01_vwn"])
    x02_vwn = parse_expr(parameters["x02_vwn"])
    x03_vwn = parse_expr(parameters["x03_vwn"])
    
    zeta = rho_s/rho_c
    rs = (3/(4*pi*rho_c))**(1/3)
    
    Q = sqrt(4*c - b**2)
    f1 = 2*b/Q
    f2 = b*x0/(x0**2 + b*x0 + c)
    f3 = 2*(2*x0 + b)/Q
    
    fpp = 4/(9*(2**(1/3) - 1 ))
    fx = rs + b*sqrt(rs) + c
    
    f_aux = A*(ln(rs/fx) + (f1 - f2*f3)*atan(Q/(2*sqrt(rs) + b)) - f2*ln((sqrt(rs) - x0)**2/fx))
    
    DMC = f_aux.subs({A:A2_vwn, b:B2_vwn, c:c2_vwn, x0:x02_vwn}) - f_aux.subs({A:A1_vwn, b:B1_vwn, c:c1_vwn, x0:x01_vwn})
    
    f_zeta = ((1 + zeta)**(4/3) + (1 - zeta)**(4/3) - 2)/(2**(4/3) - 2)
    
    eps = f_aux.subs({A:A1_vwn, b:B1_vwn, c:c1_vwn, x0:x01_vwn}) 
    eps += f_aux.subs({A:A3_vwn, b:B3_vwn, c:c3_vwn, x0:x03_vwn})*f_zeta*(1 - zeta**4)/fpp
    eps += DMC*f_zeta*zeta**4
    
    return eps
    
    
def PBEc(parameters):
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    gamma_cc = symbols('gamma_cc', real=True, nonnegative=True)
    
    # Parameters
    gamma = parse_expr(parameters["gamma_PBEc"])#gamma = (1 - ln(2))/(pi**2)
    beta = parse_expr(parameters["beta_PBEc"])#beta = 0.066725
    #beta = 0.06672455060314922 # libxc value
    
    # Functional
    e = 1
    a0 = 1
    kF = (3 * rho_c * pi**2)**(1/3)
    ks = sqrt( 4 * kF/(pi*a0) )
    zeta = rho_s/rho_c
    phi = ( (1 + zeta)**(2/3) + (1 - zeta)**(2/3) )/2
    t = sqrt(gamma_cc)/( 2 * phi * ks * rho_c )
    t2 = t**2
    ecunif = PW92c(parameters)
    A = beta/gamma*(exp(-ecunif/(gamma*phi**3*e**2/a0))-1)**(-1)
    H = (e**2/a0)*gamma*phi**3*ln(1 + beta/gamma*t**2*(( 1 + A * t**2 )/( 1 + A * t**2 + A**2 * t**4 )))
    epsPBE = ecunif + H
    
    return epsPBE
    
    
def TPSSc_case_1(parameters):
    # case 1, PBEc_alpha and PBEc_beta
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    gamma_cc, gamma_ss = symbols('gamma_cc gamma_ss', real=True, nonnegative=True)
    gamma_cs = symbols('gamma_cs', real=True)
    tau_c, tau_s = symbols('tau_c tau_s', real=True)
    
    zeta = rho_s/rho_c
    tauW = gamma_cc / (8*rho_c)
    rho_a = 1/2*(rho_c + rho_s)
    rho_b = 1/2*(rho_c - rho_s)
    
    c1 = parse_expr(parameters["c1_TPSSc"])#c1 = 0.53
    c2 = parse_expr(parameters["c2_TPSSc"])#c2 = 0.87
    c3 = parse_expr(parameters["c3_TPSSc"])#c3 = 0.5
    c4 = parse_expr(parameters["c4_TPSSc"])#c4 = 2.26
    d = parse_expr(parameters["d_TPSSc"])#d = 2.8
    
    Czero = c1 + c2*zeta**2 + c3*zeta**4 + c4*zeta**6
    abs_grad_zeta = sqrt(1/rho_c**4 * ( gamma_ss*rho_c**2 + rho_s**2*gamma_cc - 2*gamma_cs*rho_c*rho_s ) )
    eta = abs_grad_zeta / (2*(3*pi**2*rho_c)**(1/3))
    C = Czero / ( 1 + eta**2 * ( (1+zeta)**(-4/3) + (1-zeta)**(-4/3) )/2 )**4
    epsbarc_alpha = spec_func.PBEc_alpha_replaced(parameters)
    epsbarc_beta = spec_func.PBEc_beta_replaced(parameters)
    epsrevPKZB = PBEc(parameters) * ( 1 + C * (tauW/tau_c)**2 ) - (1 + C )*(tauW/tau_c)**2*(rho_a*epsbarc_alpha + rho_b*epsbarc_beta)/rho_c
    
    eps = epsrevPKZB * ( 1 + d*epsrevPKZB*(tauW/tau_c)**3 )

    return eps
    
    
def TPSSc_case_2(parameters):
    # case 2, PBEc_alpha and PBEc
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    gamma_cc, gamma_ss = symbols('gamma_cc gamma_ss', real=True, nonnegative=True)
    gamma_cs = symbols('gamma_cs', real=True)
    tau_c, tau_s = symbols('tau_c tau_s', real=True)
    
    zeta = rho_s/rho_c
    tauW = gamma_cc / (8*rho_c)
    rho_a = 1/2*(rho_c + rho_s)
    rho_b = 1/2*(rho_c - rho_s)
    
    c1 = parse_expr(parameters["c1_TPSSc"])#c1 = 0.53
    c2 = parse_expr(parameters["c2_TPSSc"])#c2 = 0.87
    c3 = parse_expr(parameters["c3_TPSSc"])#c3 = 0.5
    c4 = parse_expr(parameters["c4_TPSSc"])#c4 = 2.26
    d = parse_expr(parameters["d_TPSSc"])#d = 2.8
    
    Czero = c1 + c2*zeta**2 + c3*zeta**4 + c4*zeta**6
    abs_grad_zeta = sqrt(1/rho_c**4 * ( gamma_ss*rho_c**2 + rho_s**2*gamma_cc - 2*gamma_cs*rho_c*rho_s ) )
    eta = abs_grad_zeta / (2*(3*pi**2*rho_c)**(1/3))
    C = Czero / ( 1 + eta**2 * ( (1+zeta)**(-4/3) + (1-zeta)**(-4/3) )/2 )**4
    epsbarc_alpha = spec_func.PBEc_alpha_replaced(parameters)
    epsbarc_beta = PBEc(parameters)
    epsrevPKZB = PBEc(parameters) * ( 1 + C * (tauW/tau_c)**2 ) - (1 + C )*(tauW/tau_c)**2*(rho_a*epsbarc_alpha + rho_b*epsbarc_beta)/rho_c
    
    eps = epsrevPKZB * ( 1 + d*epsrevPKZB*(tauW/tau_c)**3 )
    
    return eps
    
    
def TPSSc_case_3(parameters):
    # case 3, PBEc and PBEc_beta
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    gamma_cc, gamma_ss = symbols('gamma_cc gamma_ss', real=True, nonnegative=True)
    gamma_cs = symbols('gamma_cs', real=True)
    tau_c, tau_s = symbols('tau_c tau_s', real=True)
    
    zeta = rho_s/rho_c
    tauW = gamma_cc / (8*rho_c)
    rho_a = 1/2*(rho_c + rho_s)
    rho_b = 1/2*(rho_c - rho_s)
    
    c1 = parse_expr(parameters["c1_TPSSc"])#c1 = 0.53
    c2 = parse_expr(parameters["c2_TPSSc"])#c2 = 0.87
    c3 = parse_expr(parameters["c3_TPSSc"])#c3 = 0.5
    c4 = parse_expr(parameters["c4_TPSSc"])#c4 = 2.26
    d = parse_expr(parameters["d_TPSSc"])#d = 2.8
    
    Czero = c1 + c2*zeta**2 + c3*zeta**4 + c4*zeta**6
    abs_grad_zeta = sqrt(1/rho_c**4 * ( gamma_ss*rho_c**2 + rho_s**2*gamma_cc - 2*gamma_cs*rho_c*rho_s ) )
    eta = abs_grad_zeta / (2*(3*pi**2*rho_c)**(1/3))
    C = Czero / ( 1 + eta**2 * ( (1+zeta)**(-4/3) + (1-zeta)**(-4/3) )/2 )**4
    epsbarc_alpha = PBEc(parameters)
    epsbarc_beta = spec_func.PBEc_beta_replaced(parameters)
    epsrevPKZB = PBEc(parameters) * ( 1 + C * (tauW/tau_c)**2 ) - (1 + C )*(tauW/tau_c)**2*(rho_a*epsbarc_alpha + rho_b*epsbarc_beta)/rho_c
    
    eps = epsrevPKZB * ( 1 + d*epsrevPKZB*(tauW/tau_c)**3 )
    
    return eps
    
    
def TPSSc_case_4(parameters):
    # case 4, PBEc and PBEc
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    gamma_cc, gamma_ss = symbols('gamma_cc gamma_ss', real=True, nonnegative=True)
    gamma_cs = symbols('gamma_cs', real=True)
    tau_c, tau_s = symbols('tau_c tau_s', real=True)
    
    zeta = rho_s/rho_c
    tauW = gamma_cc / (8*rho_c)
    rho_a = 1/2*(rho_c + rho_s)
    rho_b = 1/2*(rho_c - rho_s)
    
    c1 = parse_expr(parameters["c1_TPSSc"])#c1 = 0.53
    c2 = parse_expr(parameters["c2_TPSSc"])#c2 = 0.87
    c3 = parse_expr(parameters["c3_TPSSc"])#c3 = 0.5
    c4 = parse_expr(parameters["c4_TPSSc"])#c4 = 2.26
    d = parse_expr(parameters["d_TPSSc"])#d = 2.8
    
    Czero = c1 + c2*zeta**2 + c3*zeta**4 + c4*zeta**6
    abs_grad_zeta = sqrt(1/rho_c**4 * ( gamma_ss*rho_c**2 + rho_s**2*gamma_cc - 2*gamma_cs*rho_c*rho_s ) )
    eta = abs_grad_zeta / (2*(3*pi**2*rho_c)**(1/3))
    C = Czero / ( 1 + eta**2 * ( (1+zeta)**(-4/3) + (1-zeta)**(-4/3) )/2 )**4
    epsbarc_alpha = PBEc(parameters)
    epsbarc_beta = PBEc(parameters)
    epsrevPKZB = PBEc(parameters) * ( 1 + C * (tauW/tau_c)**2 ) - (1 + C )*(tauW/tau_c)**2*(rho_a*epsbarc_alpha + rho_b*epsbarc_beta)/rho_c
    
    eps = epsrevPKZB * ( 1 + d*epsrevPKZB*(tauW/tau_c)**3 )
    
    return eps
    
    
def LDAx(parameters):
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    
    e = 1
    rho = 2*rho_a
    kF = (3 * rho * pi**2)**(1/3)
    
    eps = -3*e**2 * kF / (4*pi)
    
    return eps
    
    
def PBEx(parameters):
    rho_a, gamma_aa = symbols('rho_a gamma_aa', real=True, nonnegative=True)
    
    rho = 2*rho_a
    gamma = 4*gamma_aa
    bPBE = parse_expr(parameters["b_PBEx"])#bPBE = 0.21951
    kappa = parse_expr(parameters["kappa_PBEx"])#kappa = 0.804
    kF = (3 * rho * pi**2)**(1/3)
    s = sqrt(gamma) / ( 2*kF*rho )
    
    Fx = 1 + kappa - kappa / ( 1 + bPBE*s**2 / kappa )
    eps = LDAx(parameters) * Fx
    
    return eps
    
    
def TPSSx(parameters):
    rho_a, gamma_aa, tau_a = symbols('rho_a gamma_aa tau_a', real=True, nonnegative=True)
    tau_a = symbols('tau_a', real=True)
    
    # Parameters
    kappa = parse_expr(parameters["kappa_PBEx"])#kkappa = 0.804
    c = parse_expr(parameters["c_TPSSx"])#c = 1.59096
    e = parse_expr(parameters["e_TPSSx"])#e = 1.537
    b = parse_expr(parameters["b_TPSSx"])#b = 0.40
    mupar = parse_expr(parameters["b_PBEx"])#mupar = 0.21951
    
    # By spin-splitting relation, rho = 2*rho_x
    rho = 2*rho_a
    # By spin-splitting relation, grad(rho) = 2*grad(rho_a), thus:
    #   gamma = 2*grad(rho_a) . 2*grad(rho_a) = 4*gamma_aa
    gamma = 4*gamma_aa
    # Since spin-splitting relation, tau = 2*tau_x
    #  and sum_sigma(tau_sigma) = tau_c = tau
    tau = 2*tau_a 
    
    p = gamma / ( 4*(3*pi**2)**(2/3) * rho**(8/3) )
    tauW = 1/8 * gamma / rho
    z = tauW / tau
    tauunif = 3/10 * (3*pi**2)**(2/3)*rho**(5/3)
    alpha = 5*p/3 * (1/z - 1)
    qbarb = 9/20 * (alpha - 1)/sqrt( 1 + b*alpha*(alpha-1) ) + 2*p/3
    x = (10/81 + c * z**2/(1 + z**2)**2)*p 
    x = x + 146/2025*qbarb**2
    x = x - 73/405*qbarb * sqrt( 1/2*(3/5 * z)**2 + 1/2*p**2 )
    x = x + 1/kappa * (10/81)**2 * p**2
    x = x + 2*sqrt(e) * 10/81 * (3/5*z)**2
    x = x + e*mupar*p**3
    x = x/(1+sqrt(e)*p)**2
    Fx = 1 + kappa - kappa / ( 1 + x/kappa )
    
    exunif = -3/(4*pi)*(3*pi**2*rho)**(1/3)
    
    # Ex = 0.5*Ex(2*rho_a) + 0.5*Ex(2*rho_b)
    eps = exunif * Fx
    
    return eps
    