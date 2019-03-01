from sympy import symbols, sqrt, exp, ln, erf, pi, log, erfc, Ei
from sympy.parsing.sympy_parser import parse_expr
import print_functional_to_DALTON as dalprint
import functionals as func


def PW92c_mu(parameters):
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    # 10.1103/PhysRevB.73.155111
    n_arg, rs_arg, mu = symbols('n_arg rs_arg mu', real=True)
    zeta = rho_s/rho_c
    rs = (3/(4*pi*rho_c))**(1/3)
    epsc = func.PW92c(parameters)
    
    phin = ( (1+zeta)**(n_arg/3) + (1-zeta)**(n_arg/3) )/2
    phi2 = phin.subs({n_arg: 2})
    phi4 = phin.subs({n_arg: 4})
    phi8 = phin.subs({n_arg: 8})
    alpha = parse_expr(parameters["alpha_srPW92c"])#alpha = (4/(9 * pi))**(1/3)
    x = mu*sqrt(rs) / (phi2)
    
    # g; from; 10.1103/PhysRevB.64.155102
    #dd = 0.7524 # Paper value
    #B = 0.73172 - dd # Paper value
    #C = 0.08193 # Paper value
    #D = -0.01277 # Paper value
    #E = 0.001859 # Paper value
    dd = parse_expr(parameters["dd_srPW92c"])#dd = 0.752411 # Dalton value
    B = parse_expr(parameters["B_srPW92c"])#B = 0.7317 - dd # Dalton value
    C = parse_expr(parameters["C_srPW92c"])#C = 0.0819306 # Dalton value
    D = parse_expr(parameters["D_srPW92c"])#D = -0.0127713 # Dalton value
    E = parse_expr(parameters["E_srPW92c"])#E = 0.00185898 # Dalton value
    #(1 - B*rs + C*rs**2 + D*rs**3 +E*rs**4)*exp(-dd*rs), there might be an argument to divide by two?
    g = (1 - B*rs + C*rs**2 + D*rs**3 +E*rs**4)*exp(-dd*rs)/2
    gc = g - 1/2
    
    a = parse_expr(parameters["a_srPW92c"])#a = 5.84605
    c = parse_expr(parameters["c_srPW92c"])#c = 3.91744
    d = parse_expr(parameters["d_srPW92c"])#d =  3.44851
    b = parse_expr(parameters["b_srPW92c"])#b = d - 3*pi*alpha/(4*ln(2)-4)
    Q = (2*ln(2) - 2)/(pi**2)*ln( (1+a*x+b*x**2+c*x**3)/(1+a*x+d*x**2) )
    
    b0 = parse_expr(parameters["b0_srPW92c"])*rs #b0 = 0.784949*rs
    bpp = parse_expr(parameters["bpp_srPW92c"])#bpp = 0.4319
    app = parse_expr(parameters["app_srPW92c"])#app = bpp - 0.454555 #-0.02267
    cpp = parse_expr(parameters["cpp_srPW92c"])#cpp = 0.04
    gpp = 2**(5/3)/(5*alpha**2 * rs_arg**2) * ( 1+app*rs_arg ) / ( 1+bpp*rs_arg + cpp*rs_arg**2 )
    D2 = exp(-0.547*rs)/(rs**2) * ( -0.388*rs + 0.676*rs**2 )
    D3 = exp(-0.31*rs)/(rs**3)*(-4.95*rs + rs**2)
    
    c4 = ((1+zeta)/2)**2*gpp.subs({rs_arg: rs*(2/(1+zeta))**(1/3)})
    c4 = c4 + ((1-zeta)/2)**2*gpp.subs({rs_arg: rs*(2/(1-zeta))**(1/3)})
    c4 = c4 + (1 - zeta**2) * D2 - phi8/(5*alpha**2*rs**2)
    c5 = ((1+zeta)/2)**2*gpp.subs({rs_arg: rs*(2/(1+zeta))**(1/3)})
    c5 = c5 + ((1-zeta)/2)**2*gpp.subs({rs_arg: rs*(2/(1-zeta))**(1/3)})
    c5 = c5 + (1 - zeta**2) * D3
    
    C2 = -( 3*(1-zeta**2)*gc ) / ( 8*rs**3 )
    C3 = -(1-zeta**2)*g/(sqrt(2*pi)*rs**3)
    C4 = -9*c4/(64*rs**3)
    C5 = -9*c5/(40*sqrt(2*pi)*rs**3)
    
    a1 = 4*b0**6*C3 + b0**8*C5
    a2 = 4*b0**6*C2 + b0**8*C4 + 6*b0**4*epsc
    a3 = b0**8*C3
    a4 = b0**8*C2 + 4*b0**6*epsc
    a5 = b0**8*epsc
    
    epsc_LR = phi2**3 * Q + a1*mu**3 + a2*mu**4 + a3*mu**5 + a4*mu**6 +a5*mu**8
    epsc_LR = epsc_LR / (( 1 + b0**2 * mu**2 )**4)
    eps = epsc - epsc_LR
    
    return eps
    
    
def VWN5c_mu(parameters):
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    # 10.1103/PhysRevB.73.155111
    n_arg, rs_arg, mu = symbols('n_arg rs_arg mu', real=True)
    zeta = rho_s/rho_c
    rs = (3/(4*pi*rho_c))**(1/3)
    epsc = func.VWN5c(parameters)
    
    phin = ( (1+zeta)**(n_arg/3) + (1-zeta)**(n_arg/3) )/2
    phi2 = phin.subs({n_arg: 2})
    phi4 = phin.subs({n_arg: 4})
    phi8 = phin.subs({n_arg: 8})
    alpha = parse_expr(parameters["alpha_srPW92c"])#alpha = (4/(9 * pi))**(1/3)
    x = mu*sqrt(rs) / (phi2)
    
    # g; from; 10.1103/PhysRevB.64.155102
    #dd = 0.7524 # Paper value
    #B = 0.73172 - dd # Paper value
    #C = 0.08193 # Paper value
    #D = -0.01277 # Paper value
    #E = 0.001859 # Paper value
    dd = parse_expr(parameters["dd_srPW92c"])#dd = 0.752411 # Dalton value
    B = parse_expr(parameters["B_srPW92c"])#B = 0.7317 - dd # Dalton value
    C = parse_expr(parameters["C_srPW92c"])#C = 0.0819306 # Dalton value
    D = parse_expr(parameters["D_srPW92c"])#D = -0.0127713 # Dalton value
    E = parse_expr(parameters["E_srPW92c"])#E = 0.00185898 # Dalton value
    #(1 - B*rs + C*rs**2 + D*rs**3 +E*rs**4)*exp(-dd*rs), there might be an argument to divide by two?
    g = (1 - B*rs + C*rs**2 + D*rs**3 +E*rs**4)*exp(-dd*rs)/2
    gc = g - 1/2
    
    a = parse_expr(parameters["a_srPW92c"])#a = 5.84605
    c = parse_expr(parameters["c_srPW92c"])#c = 3.91744
    d = parse_expr(parameters["d_srPW92c"])#d =  3.44851
    b = parse_expr(parameters["b_srPW92c"])#b = d - 3*pi*alpha/(4*ln(2)-4)
    Q = (2*ln(2) - 2)/(pi**2)*ln( (1+a*x+b*x**2+c*x**3)/(1+a*x+d*x**2) )
    
    b0 = parse_expr(parameters["b0_srPW92c"])*rs #b0 = 0.784949*rs
    bpp = parse_expr(parameters["bpp_srPW92c"])#bpp = 0.4319
    app = parse_expr(parameters["app_srPW92c"])#app = bpp - 0.454555 #-0.02267
    cpp = parse_expr(parameters["cpp_srPW92c"])#cpp = 0.04
    gpp = 2**(5/3)/(5*alpha**2 * rs_arg**2) * ( 1+app*rs_arg ) / ( 1+bpp*rs_arg + cpp*rs_arg**2 )
    D2 = exp(-0.547*rs)/(rs**2) * ( -0.388*rs + 0.676*rs**2 )
    D3 = exp(-0.31*rs)/(rs**3)*(-4.95*rs + rs**2)
    
    c4 = ((1+zeta)/2)**2*gpp.subs({rs_arg: rs*(2/(1+zeta))**(1/3)})
    c4 = c4 + ((1-zeta)/2)**2*gpp.subs({rs_arg: rs*(2/(1-zeta))**(1/3)})
    c4 = c4 + (1 - zeta**2) * D2 - phi8/(5*alpha**2*rs**2)
    c5 = ((1+zeta)/2)**2*gpp.subs({rs_arg: rs*(2/(1+zeta))**(1/3)})
    c5 = c5 + ((1-zeta)/2)**2*gpp.subs({rs_arg: rs*(2/(1-zeta))**(1/3)})
    c5 = c5 + (1 - zeta**2) * D3
    
    C2 = -( 3*(1-zeta**2)*gc ) / ( 8*rs**3 )
    C3 = -(1-zeta**2)*g/(sqrt(2*pi)*rs**3)
    C4 = -9*c4/(64*rs**3)
    C5 = -9*c5/(40*sqrt(2*pi)*rs**3)
    
    a1 = 4*b0**6*C3 + b0**8*C5
    a2 = 4*b0**6*C2 + b0**8*C4 + 6*b0**4*epsc
    a3 = b0**8*C3
    a4 = b0**8*C2 + 4*b0**6*epsc
    a5 = b0**8*epsc
    
    epsc_LR = phi2**3 * Q + a1*mu**3 + a2*mu**4 + a3*mu**5 + a4*mu**6 +a5*mu**8
    epsc_LR = epsc_LR / (( 1 + b0**2 * mu**2 )**4)
    eps = epsc - epsc_LR
    
    return eps
    
    
def PBEc_mu(parameters):
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    gamma_cc = symbols('gamma_cc', real=True, nonnegative=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    betaPBE = parse_expr(parameters["beta_PBEc"])#betaPBE = 0.066725
    gamma = parse_expr(parameters["gamma_PBEc"])#gamma = 0.031091
    alpha_c = parse_expr(parameters["alphac_srPBEc"])#alpha_c = 2.78
    
    zeta = rho_s/rho_c
    phi = ( (1 + zeta)**(2/3) + (1 - zeta)**(2/3) )/2
    
    kF = (3*pi**2*rho_c)**(1/3)
    ks = sqrt(4*kF/pi)
    t = sqrt(gamma_cc)/(2*phi*ks*rho_c)
    
    beta = betaPBE*(PW92c_mu(parameters)/func.PW92c(parameters))**alpha_c
    A = beta / ( gamma * ( exp(-PW92c_mu(parameters)/(gamma*phi**3)) -1 ) )
    H = gamma*phi**3 * ln( 1 + beta*t**2/gamma * ( (1+A*t**2)/(1+A*t**2+A**2*t**4) ) )
    
    eps = PW92c_mu(parameters) + H
    
    return eps

    
def TPSSc_mu_case_1(parameters):
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    gamma_cc, gamma_ss = symbols('gamma_cc gamma_ss', real=True, nonnegative=True)
    gamma_cs = symbols('gamma_cs', real=True)
    tau_c, tau_s = symbols('tau_c tau_s', real=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    etac = parse_expr(parameters["etac_srTPSSc"])#etac = 2.9
    kF = (3*pi**2*rho_c)**(1/3)
    mubar = mu / ( 2*kF )
    
    eps = PBEc_mu(parameters) + (func.TPSSc_case_1(parameters) - func.PBEc(parameters))*exp(-etac*mubar)
    
    return eps
    
    
def TPSSc_mu_case_2(parameters):
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    gamma_cc, gamma_ss = symbols('gamma_cc gamma_ss', real=True, nonnegative=True)
    gamma_cs = symbols('gamma_cs', real=True)
    tau_c, tau_s = symbols('tau_c tau_s', real=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    etac = parse_expr(parameters["etac_srTPSSc"])#etac = 2.9
    kF = (3*pi**2*rho_c)**(1/3)
    mubar = mu / ( 2*kF )
    
    eps = PBEc_mu(parameters) + (func.TPSSc_case_2(parameters) - func.PBEc(parameters))*exp(-etac*mubar)
    
    return eps
    
    
def TPSSc_mu_case_3(parameters):
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    gamma_cc, gamma_ss = symbols('gamma_cc gamma_ss', real=True, nonnegative=True)
    gamma_cs = symbols('gamma_cs', real=True)
    tau_c, tau_s = symbols('tau_c tau_s', real=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    etac = parse_expr(parameters["etac_srTPSSc"])#etac = 2.9
    kF = (3*pi**2*rho_c)**(1/3)
    mubar = mu / ( 2*kF )
    
    eps = PBEc_mu(parameters) + (func.TPSSc_case_3(parameters) - func.PBEc(parameters))*exp(-etac*mubar)
    
    return eps
    
    
def TPSSc_mu_case_4(parameters):
    rho_c = symbols('rho_c', real=True, nonnegative=True)
    rho_s = symbols('rho_s', real=True)
    gamma_cc, gamma_ss = symbols('gamma_cc gamma_ss', real=True, nonnegative=True)
    gamma_cs = symbols('gamma_cs', real=True)
    tau_c, tau_s = symbols('tau_c tau_s', real=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    etac = parse_expr(parameters["etac_srTPSSc"])#etac = 2.9
    kF = (3*pi**2*rho_c)**(1/3)
    mubar = mu / ( 2*kF )
    
    eps = PBEc_mu(parameters) + (func.TPSSc_case_4(parameters) - func.PBEc(parameters))*exp(-etac*mubar)
    
    return eps
    

def LDAx_mu_case_1(parameters):
    # using LDAx from Dalton
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    rho = 2*rho_a
    AKF = rho**(1/3)*((3*pi**2)**(1/3))
    A = mu/(2*AKF)
    
    # Case 1 - A .lt. 10**-9
    eps = -3/8 * (24*rho/pi)**(1/3)
    
    return eps
    
    
def LDAx_mu_case_2(parameters):
    # using LDAx from Dalton
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    rho = 2*rho_a
    AKF = rho**(1/3)*((3*pi**2)**(1/3))
    A = mu/(2*AKF)
    
    # Case 2 - A .le. 100    
    eps = - ((24*rho/pi)**(1/3))*(3/8 - A*(sqrt(pi)*erf(0.5/A) +  (2*A - 4*A**3)*exp(-0.25/A**2) - 3*A + 4*A**3))
    
    return eps
    
    
def LDAx_mu_case_3(parameters):
    # using LDAx from Dalton
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    rho = 2*rho_a
    AKF = rho**(1/3)*((3*pi**2)**(1/3))
    A = mu/(2*AKF)
    
    # Case 3 - A .lt. 10**9
    eps = - ((24*rho/pi)**(1/3)) * 1 / (96*A**2)
    
    return eps
    
    
def PBEx_mu_case_1(parameters):
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    gamma_aa = symbols('gamma_aa', real=True, nonnegative=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    kappa = parse_expr(parameters["kappa_PBEx"])#kappa = 0.804
    bPBE = parse_expr(parameters["b_PBEx"])#bPBE = 0.21951
    alphax = parse_expr(parameters["alphax_srPBEx"])#alphax = 19
    rho = 2*rho_a
    gamma = 4*gamma_aa
    kF = (3*pi**2*rho)**(1/3)
    a = mu/(2*kF)
    mubar = mu / ( 2*kF )
    fak=2.540118935556*exp(-alphax*a**2)
    s = sqrt(gamma) / ( 2*kF*rho )
    
    # Case 1 - A .lt. 10**-9
    # can only have berf case 1
    # using berf from Dalton
    berf = ((-7+72*a**2)/(27*(-3-24*a**2+32*a**4+8*sqrt(pi)*a)))*fak
    Fx = 1 + kappa - kappa / ( 1 + berf * s**2 / kappa)
    eps = LDAx_mu_case_1(parameters) * Fx
    
    return eps
    
    
def PBEx_mu_case_2_1(parameters):
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    gamma_aa = symbols('gamma_aa', real=True, nonnegative=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    kappa = parse_expr(parameters["kappa_PBEx"])#kappa = 0.804
    bPBE = parse_expr(parameters["b_PBEx"])#bPBE = 0.21951
    alphax = parse_expr(parameters["alphax_srPBEx"])#alphax = 19
    rho = 2*rho_a
    gamma = 4*gamma_aa
    kF = (3*pi**2*rho)**(1/3)
    a = mu/(2*kF)
    mubar = mu / ( 2*kF )
    fak=2.540118935556*exp(-alphax*a**2)
    s = sqrt(gamma) / ( 2*kF*rho )
    
    # Case 2 - A .le. 100
    # can have berf case 1, 2 or 3
    # using berf from Dalton
    # Case 2.1 - A .lt. 0.075
    berf = ((-7+72*a**2)/(27*(-3-24*a**2+32*a**4+8*sqrt(pi)*a)))*fak
    Fx = 1 + kappa - kappa / ( 1 + berf * s**2 / kappa)
    eps = LDAx_mu_case_2(parameters) * Fx
    
    return eps
    
    
def PBEx_mu_case_2_2(parameters):
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    gamma_aa = symbols('gamma_aa', real=True, nonnegative=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    kappa = parse_expr(parameters["kappa_PBEx"])#kappa = 0.804
    bPBE = parse_expr(parameters["b_PBEx"])#bPBE = 0.21951
    alphax = parse_expr(parameters["alphax_srPBEx"])#alphax = 19
    rho = 2*rho_a
    gamma = 4*gamma_aa
    kF = (3*pi**2*rho)**(1/3)
    a = mu/(2*kF)
    mubar = mu / ( 2*kF )
    fak=2.540118935556*exp(-alphax*a**2)
    s = sqrt(gamma) / ( 2*kF*rho )
    
    # Case 2 - A .le. 100
    # can have berf case 1, 2 or 3
    # using berf from Dalton
    # Case 2.2 - A .gt. 50
    berf = (1/(72*a**2)-1/(17280*a**4) - 23/(358400*a**6))*fak
    Fx = 1 + kappa - kappa / ( 1 + berf * s**2 / kappa)
    eps = LDAx_mu_case_2(parameters) * Fx
    
    return eps
    
    
def PBEx_mu_case_2_3(parameters):
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    gamma_aa = symbols('gamma_aa', real=True, nonnegative=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    kappa = parse_expr(parameters["kappa_PBEx"])#kappa = 0.804
    bPBE = parse_expr(parameters["b_PBEx"])#bPBE = 0.21951
    alphax = parse_expr(parameters["alphax_srPBEx"])#alphax = 19
    rho = 2*rho_a
    gamma = 4*gamma_aa
    kF = (3*pi**2*rho)**(1/3)
    a = mu/(2*kF)
    mubar = mu / ( 2*kF )
    fak=2.540118935556*exp(-alphax*a**2)
    s = sqrt(gamma) / ( 2*kF*rho )
    
    # Case 2 - A .le. 100
    # can have berf case 1, 2 or 3
    # using berf from Dalton
    # Case 2.3 - A between 0.075 and 50
    fak2 = exp(0.25/a**2)
    berf = (1.851851851851851851851852*10**-2*(-1 + 1.44*10**2*a**4*(-1 + fak2) - 2*a**2*(1.1*10 + 7*fak2 ))) / (a**2*(3.2*10*a**4*(-1 + fak2) - 3*fak2 + 1.417963080724412821838534*10*a*erf(5*10**-1/a)*fak2 - 8*a**2*(-2 + 3*fak2)))
    berf *= fak
    Fx = 1 + kappa - kappa / ( 1 + berf * s**2 / kappa)
    eps = LDAx_mu_case_2(parameters) * Fx
    
    return eps
    
    
def PBEx_mu_case_3(parameters):
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    gamma_aa = symbols('gamma_aa', real=True, nonnegative=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    kappa = parse_expr(parameters["kappa_PBEx"])#kappa = 0.804
    bPBE = parse_expr(parameters["b_PBEx"])#bPBE = 0.21951
    alphax = parse_expr(parameters["alphax_srPBEx"])#alphax = 19
    rho = 2*rho_a
    gamma = 4*gamma_aa
    kF = (3*pi**2*rho)**(1/3)
    a = mu/(2*kF)
    mubar = mu / ( 2*kF )
    fak=2.540118935556*exp(-alphax*a**2)
    s = sqrt(gamma) / ( 2*kF*rho )
    
    # Case 3 - A .lt. 10**9
    # can only have berf case 2
    # using berf from Dalton
    berf = (1/(72*a**2)-1/(17280*a**4) - 23/(358400*a**6))*fak
    Fx = 1 + kappa - kappa / ( 1 + berf * s**2 / kappa)
    eps = LDAx_mu_case_3(parameters) * Fx
    
    return eps


def TPSSx_mu_case_1(parameters):
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    gamma_aa = symbols('gamma_aa', real=True, nonnegative=True)
    tau_a = symbols('tau_a', real=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    rho = 2*rho_a
    etax = parse_expr(parameters["etax_srTPSSx"])#etax = 15
    kF = (3*pi**2*rho)**(1/3)
    mubar = mu / ( 2*kF )
    
    eps = PBEx_mu_case_1(parameters) + (func.TPSSx(parameters) - func.PBEx(parameters))*exp(-etax*mubar)

    return eps
    
    
def TPSSx_mu_case_2_1(parameters):
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    gamma_aa = symbols('gamma_aa', real=True, nonnegative=True)
    tau_a = symbols('tau_a', real=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    rho = 2*rho_a
    etax = parse_expr(parameters["etax_srTPSSx"])#etax = 15
    kF = (3*pi**2*rho)**(1/3)
    mubar = mu / ( 2*kF )
    
    eps = PBEx_mu_case_2_1(parameters) + (func.TPSSx(parameters) - func.PBEx(parameters))*exp(-etax*mubar)
    
    return eps
    
    
def TPSSx_mu_case_2_2(parameters):
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    gamma_aa = symbols('gamma_aa', real=True, nonnegative=True)
    tau_a = symbols('tau_a', real=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    rho = 2*rho_a
    etax = parse_expr(parameters["etax_srTPSSx"])#etax = 15
    kF = (3*pi**2*rho)**(1/3)
    mubar = mu / ( 2*kF )
    
    eps = PBEx_mu_case_2_2(parameters) + (func.TPSSx(parameters) - func.PBEx(parameters))*exp(-etax*mubar)
    
    return eps
    
    
def TPSSx_mu_case_2_3(parameters):
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    gamma_aa = symbols('gamma_aa', real=True, nonnegative=True)
    tau_a = symbols('tau_a', real=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    rho = 2*rho_a
    etax = parse_expr(parameters["etax_srTPSSx"])#etax = 15
    kF = (3*pi**2*rho)**(1/3)
    mubar = mu / ( 2*kF )
    
    eps = PBEx_mu_case_2_3(parameters) + (func.TPSSx(parameters) - func.PBEx(parameters))*exp(-etax*mubar)
    
    return eps
    
    
def TPSSx_mu_case_3(parameters):
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    gamma_aa = symbols('gamma_aa', real=True, nonnegative=True)
    tau_a = symbols('tau_a', real=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    rho = 2*rho_a
    etax = parse_expr(parameters["etax_srTPSSx"])#etax = 15
    kF = (3*pi**2*rho)**(1/3)
    mubar = mu / ( 2*kF )
    
    eps = PBEx_mu_case_3(parameters) + (func.TPSSx(parameters) - func.PBEx(parameters))*exp(-etax*mubar)
    
    return eps


def wPBEx(parameters):
    rho_a = symbols('rho_a', real=True, nonnegative=True)
    gamma_aa = symbols('gamma_aa', real=True, nonnegative=True)
    mu = symbols('mu', real=True, nonnegative=True)
    
    Abar = 0.757211
    D = 0.609650
    B = 9/4*(Abar*D - 1/2*Abar**2) - 1/2 # -0.106364
    C = 1/10 - 9/8*Abar*(D**2 - 1/3*Abar**2) + B*D # -0.118649
    #E = -2/5*C*D- 4/15*B*D**2 - 4/5*sqrt(pi)*D**(7/2) + 6/5*sqrt(Abar)*D**3*(2*sqrt(D) - sqrt(Abar)) # -0.0477963
    s0 = 2
    
    a2 = 0.0159941
    a3 = 0.0852995
    a4 = -0.160368
    a5 = 0.152645
    a6 = -0.0971263
    a7 = 0.0422061
    b1 = 5.33319
    b2 = -12.4780
    b3 = 11.0988
    b4 = -5.11013
    b5 = 1.71468
    b6 = -0.610380
    b7 = 0.307555
    b8 = -0.0770547
    b9 = 0.0334840
    
    rho = 2*rho_a
    gamma = 4*gamma_aa
    
    kF = (3*pi**2*rho)**(1/3)
    s = sqrt(gamma)/(2*kF*rho)
    
    H = a2*s**2 + a3*s**3 + a4*s**4 + a5*s**5 + a6*s**6 + a7*s**7
    H /= 1 + b1*s + b2*s**2 + b3*s**3 + b4*s**4 +b5*s**5 + b6*s**6 + b7*s**7 + b8*s**8 + b9*s**9
    
    Fbar = 1 - 1/(27*C) * s**2/( 1 + s**2/s0**2 ) - 1/(2*C)*s**2*H
    
    zeta = s**2*H
    eta = Abar + s**2*H
    Lambda = D + s**2*H
    
    # Gbar == E * Gbar
    Gbar = -2/5*C*Fbar*Lambda
    Gbar += -4/15*B*Lambda**2
    Gbar += -6/5*Abar*Lambda**3
    Gbar += -4/5*sqrt(pi)*Lambda**(7/2)
    Gbar += -12/5*Lambda**(7/2)*(sqrt(zeta) - sqrt(eta))
    
    nu = mu/kF
    chi = nu / sqrt( Lambda + nu**2 )
    
    F_wPBE = Abar - 4/9 * B / Lambda * (1 - chi)
    F_wPBE += -4/9 * C*Fbar / Lambda**2 * ( 1 - 3/2 * chi + 1/2 * chi**3 )
    # Gbar == E * Gbar
    F_wPBE += -8/9 * Gbar / Lambda**3 * ( 1 - 15/8 * chi + 5/4 * chi**3 - 3/8 * chi**5 )
    F_wPBE += 2*nu * ( sqrt( zeta + nu**2 ) - sqrt( eta + nu**2 ) )
    F_wPBE += 2*zeta * ln( ( nu + sqrt( zeta + nu**2 ) )/( nu + sqrt( Lambda + nu**2 ) ) )
    F_wPBE += -2*eta * ln( ( nu + sqrt( eta + nu**2) )/( nu + sqrt( Lambda + nu**2 ) ) )
    
    eps = func.LDAx(parameters) * F_wPBE
    
    return eps
    
    
