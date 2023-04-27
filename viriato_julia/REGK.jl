using LinearAlgebra 
using .Brackets, .constants,.Diag,.Aux # Ought to define these as "modules" using module __ at beginning, then choosing what data to export



# Define/Import Stuff, Allocate arrays here not sure if optimal?

# phik,nek,akpar,uekpar,gk,dummy_real,bracket_akpar_uekpar

g_inc = true # would be in constants file!
# Begin setup of necessary values, dx,dys, hypercoeffs, timestep

# Take derivaties necesary to calculate hypercoeffs, timestep for first iteration
dxphi,dyphi = convol(phik) # Better to write functions like this or as something that modifies value?
dxne, dyne = convol(nek)
dxapar,dyapar = convol(akpar)
dxuepar,dyuepar = convol(uekpar)

if g_inc 
    for ng = 1:ngtot
        dxg[:,:,ng],dyg[:,:,ng] = convol(gk[:,:,ng])
    end
end 

vex,vey,vrhosx,vrhosy = flows(dxphi,dyphi)
vxmax = maximum(max(abs(vex),abs(vrhosx)))
vymax = maximum(max(abs(vey),abs(vrhosxy)))

bx,by = bfield(dxapar,dyapar)
bxmax = maximum(abs(bx))
bymax = maxumum(abs(by))
bperp = sqrt.(bx.^2+by.^2)
bperp_max = maximum(bperp)

omega_kaw(bperp_max) # Function omega_kaw--its value is public in the function def

# Calculate CFL fraction/timestep
if g_inc
    CFL_flow= min(dx/vxmax,dy/vymax,2.0/omega_kaw,
    (1.0/rhos_de)*1.0/sqrt(ngtot*1.0)*min(dx/bxmax,dy/bymax)) # The other lines have to do with prop in z direction
else 
    CFL_flow=min(dx/vxmax,dy/vymax,dy/bymax,dx/bxmax,2.0/omega_kaw)
end

dti = CFL_frac*CFL_flow # TIMESTEP!

# Hypercoefficients!
if hyper_fixed 
    ν_g = hyper_ν_g
    ν2 = hyper_ν
    η2 = hyper_η
else
    ν_g = hyper_coef_g/dti/kperpmax^(2*hyper_order_g)
    ν2 = hyper_coef/dti/kperpmax^(2*hyper_order)
    if kperpmax^2*de^2 > 1
        η2 = hyper_coef/dti/kperpmax^(2*hyper_order-2)*de^2
    else
        res2=hyper_coef/dti/kperpmax^(2*hyper_order)
    end
end

if hyper_colls_fixed
    hyper_νei=hyper_colls
else    
    hyper_νei=hyperm_coef/dti/(Ngtot+1)^(2*hyper_morder)
end


############## TIME LOOP ###################

#First some values which must be initialized for LOOP
relative_error = 0.0
aa0 = init_aa0_fac 

p = 0
z = 0 # not necessary if 2d? Not sure if this is z as in z direction
repeat = false # Flags for certain behaviors in loop. Trying again with smaller timestep, for example
noinc = false
divergent = false
first = true # in priciple this should not be true if restarts are enabled, but lets get to that later
for t = 0:tmax 
    p = t-z # ? 

    if repeat 
        repeat = false
    else
        if divergent
            divergent = false
        else

            p_count = 0 # Number of loops through corrector step 
            if first 

                dxphi,dyphi = convol(phik)
                dxne, dyne = convol(nek)
                dxapar,dyapar = convol(akpar)
                dxuepar,dyuepar = convol(uekpar)

                if g_inc 
                    for ng = 1:ngtot
                        dxg[:,:,ng],dyg[:,:,ng] = convol(gk[:,:,ng])
                    end
                end 
            end
            first = false

            # Get the nonlinear operator values
            
            #\mathcal{N}
            fne_old, bracket_akpar_uekpar = func_ne(dxphi,dyphi,dxne,dyne,dxapar,dyapar,dxuepar,dyuepar)
            
            if g_inc
                # \mathcal{A}
                fApar_old = func_Akpar(dxapar,dyapar,dxphi,dyphi,dxne,dyne,dxuepar,dyuepar,dxg[:,:,gmin],dyg[:,:,gmin])
                # \mathcal{g2}
                fg2_old = func_g2(dxg[:,:,gmin],dyg[:,:,gmin],dxphi,dyphi,dxapar,dyapar,
                dxg[:,:,gmin+1],dyg[:,:,gmin+1],phik)
                # \mathcal{gm}
                for ng in gmin+1:ngtot-1
                    fgm_old[:,:,ng] = funcgm(ng,dxg[:,:,ng-1],dyg[:,:,ng-1],dxg[:,:,ng],dyg[:,:,ng],dxg[:,:,ng+1],dyg[:,:,ng+1],
                    dxphi,dyphi,dxapar,dyapar,akpar) 
                end
                # \mathcal{glast}
                fglast_old = func_lastg(hyper_νei,η2,dxg[:,:,ngtot-1],dyg[:,:,ngtot-1],dxg[:,:,ngtot],dyg[:,:,ngtot],
                dxphi,dyphi,dxapar,dyapar)

            else 
                funcAkpar(dxapar,dyapar,dxphi,dyphi,dxne,dyne,dxuepar,dyuepar,dummy_real,dummy_real)
            end
        end 
    end

    # Get SI operator necessary for corrector step loop
    semi_implicit_operator = calc_semi_implicit_operator(dti,bperp_max,aa0)

    # Start Predictor ("star") step
    
