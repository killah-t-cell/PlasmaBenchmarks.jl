module PlasmaBenchmarks

# TODO
# add PIC benchmarks

@info "compiling Plasma"
using Plasma
using BSON

@info "starting 1D1V model"
TD = 30000 # eV
Te = 10000 # eV

D = species.D
e = species.e

function HotCarrier(T) 
    Kb = 8.617333262145e-5
    function P(x,v)
        v_ = sqrt(sum(v .^2))
        exp(-v_/(Kb*T))
    end
end

D_D = Distribution(HotCarrier(TD), D)
G = Geometry() 

plasma = ElectrostaticPlasma([D_D], G)

sol1 = Plasma.solve(plasma, dim=1, GPU=false, inner_layers=32, time_ub=100.0, ub=10.0; strategy=QuadratureTraining()) 
sol2 = Plasma.solve(plasma, dim=1, GPU=false, inner_layers=32, time_ub=100.0, ub=10.0; strategy=StochasticTraining(200)) 

@save "1d1v_quadrature.bson" sol1
@save "1d1v_quadrature.bson" sol2

Plasma.plot(sol1)
Plasma.plot(sol2)


#########################
@info "starting 3D3V model"
TD = 15000 # eV
TT = 15000 # eV
Te = 13000 # eV

e = species.e
T = species.T
D = species.D

De = Distribution(Maxwellian(Te, e.m), e)
DT = Distribution(Maxwellian(TT, T.m), T)
DD = Distribution(Maxwellian(TD, D.m), D)
G = Geometry()

plasma = CollisionlessPlasma([De,DT,DD], G)

sol = Plasma.solve(plasma, GPU=false)

@save "3d3v.bson" sol

Plasma.plot(sol)

end

#########################
@info "starting 3D3V electrostatic model"
TD = 70000 # eV

D = species.D

DD = Distribution(HotCarrier(TD), D)
G = Geometry()

plasma = ElectrostaticPlasma([DD], G)

sol = Plasma.solve(plasma; GPU=false, time_ub=10.0)

@save "3d3v_electrostatic.bson" sol

Plasma.plot(sol)

end