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