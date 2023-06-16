include("../assets/params.jl")

function pk_pd!(du, u, p, t)
    #unpack parameters
    (gamma_1,psi,C0,D0,r,K,BW,IC50_1,Imax_1,IC50_2,gamma_2,Imax_2,xi,VD1,Cl1,k23,ka1,k32,Cl2,ka2,Vpla,Q,Vtis) = p[1:length(ode_params)]
    #unpack state variables
    C, D, AbsTMZ, PlaTMZ, CSFTMZ, AbsRG, PlaRG, TisRG= u

    #concentrations for the drug effects below
    cCSFTMZ = CSFTMZ / 140
    cPlaRG = PlaRG / (1000 * Vpla)

    #combination drug effects function
    exp1 = real(complex(cCSFTMZ)^gamma_1)/((psi*IC50_1)^gamma_1)
    exp2 = real(complex(xi*cPlaRG)^gamma_2)/((psi*IC50_1)^gamma_2)
    E = ((Imax_1*exp1)+(Imax_2*exp2)+((Imax_1+Imax_2-(Imax_1*Imax_2))*exp1*exp2))/(exp1+exp2+(exp1*exp2)+1) 

    #efect on tumor growth
    t1=-1/r*log(log(C/K)/log(C0/K))
    t2=t1+3*24
    fun = K*exp(log(C0/K)*exp(-r*t2))
    delta=(E/72)*(fun/C)

    dC = C*r*log(K/C)-delta*C
    dD = delta*C

    #drug PK
    dAbsTMZ = -ka1 * AbsTMZ
    dPlaTMZ = ka1 * AbsTMZ - Cl1 / VD1 * PlaTMZ - k23 * PlaTMZ + k32 * CSFTMZ
    dCSFTMZ = k23 * PlaTMZ - k32 * CSFTMZ

    dAbsRG = -ka2 * AbsRG
    dPlaRG = ka2*AbsRG-(Cl2/Vpla)*PlaRG+(Q/Vtis)*TisRG-(Q/Vpla)*PlaRG
    dTisRG = -(Q/Vtis)*TisRG+(Q/Vpla)*PlaRG

    #pack rhs
    du .= [dC, dD, dAbsTMZ, dPlaTMZ, dCSFTMZ, dAbsRG, dPlaRG, dTisRG]
end