# https://selalib.github.io/vlasov-maxwell.html
# see reference here for conditions

@info "starting Weibel model"

e = species.e

function Weibel(vth, k, β, α, Tr) 

    function P(x,v)
        if !(v isa Array)
            v = [v]    
        end

        if !(x isa Array)
            x = [x]    
        end

        v = sum(v .^2)
        x = sum(x .^2)

        1/(π*vth*sqrt(Tr)) * exp(-(v/Tr) / vth^2) * (1 + α*cos(k*x))
    end
end

D_e = Distribution(Weibel(0.02, 1.25, 10e-4, 10e-4, 48), e) 
G = Geometry() 

plasma = ElectrostaticPlasma([D_e], G)

sol = Plasma.solve(plasma, dim=3, GPU=false, time_ub = 50.0, ub=2*π) 

@save "weibel.bson" sol

Plasma.plot(sol)