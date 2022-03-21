# PlasmaBenchmarks

## Two stream instability

```julia
using Plasma

e = species.e

function TwoStream(vs1, vs2) 

    function P(x,v)
        if !(v isa Array)
            v = [v]    
        end

        if !(x isa Array)
            x = [x]    
        end

        v = sqrt(sum(v .^2))
        x = sqrt(sum(x .^2))

        1/2*(1/(sqrt(2*π)*vs1)*exp(-(v-vs2)^2)/(2*vs1^2) + 1/(sqrt(2*π)*vs1)*exp(-(v+vs2)^2)/(2*vs1^2))
    end
end

D_e = Distribution(TwoStream(1.6, -1.4), e) 
G = Geometry() 

plasma = ElectrostaticPlasma([D_e], G)

sol = Plasma.solve(plasma, dim=1, GPU=false, time_ub = 4.0, ub=4.0, lb=-4.0) 
```
