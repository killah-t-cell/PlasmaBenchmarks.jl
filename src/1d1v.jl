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