using LinearAlgebra, FFTW, Logging

include("constants.jl")
include("transforms.jl")
include("Brackets.jl")
include("aux.jl")
include("diag.jl")
include("functions.jl")
include("initialize.jl")
# For now, ignore the stuff that has to do with turb, anjor, 3d, antenna

# Define/Import Stuff, Allocate arrays here not sure if optimal?

# Only allocating the arrays that would be undefined. Many are returned from a function before use, for 
# example, dxapar, etc. However, dxg must be indexed the first time it is used, so must define here. Not sure if this is the best way

# LIST OF NEEDED CONSTANTS:

# nkx,nky,ngtot,gmin,ginc,dx,dy,CFL_frac,hyper_ν_g, hyper_ν, hyper_η,rhos_de,hyper_order,hyper_order_g
# hyper_coef,hyper_coef_g,kperpmax,hyperm_coef,hyper_morder, init_aa0_fac,tmax, low

function main() # This seems to be a way to reduce use of "global"

#vincent: I added the logger here. for messages that are debug messages you want to log, just do @debug "message". For messages that are of "higher importance", such as maybe logging the parameters (I am not that familiar with what)
# should actually be logged, use @info "parameter". @info is a higher importance than @debug, so if you change it to logger = Logging.SimpleLogger(log_file, Logging.Info) it will only record the @Info messages
log_file = open("mylog.log", "w")
logger = Logging.SimpleLogger(log_file, Logging.Debug)

Logging.global_logger(logger)


debugging = true # Debugging print statements


phik = Array{ComplexF64}(undef,nkx,nky)
nek = Array{ComplexF64}(undef,nkx,nky)
akpar = Array{ComplexF64}(undef,nkx,nky)
uekpar = Array{ComplexF64}(undef,nkx,nky)
nek_perturb = Array{ComplexF64}(undef,nkx,nky)

dummy_real = Array{Float64}(undef,nlx,nly)

gk = Array{ComplexF64}(undef,nkx,nky,ngtot)
dxg = Array{Float64}(undef,nlx,nly,ngtot)
dyg = Array{Float64}(undef,nlx,nly,ngtot)
fgm_old = Array{ComplexF64}(undef,nkx,nky,ngtot)

nek_star = Array{ComplexF64}(undef,nkx,nky)
phik_star = Array{ComplexF64}(undef,nkx,nky)
akpar_star = Array{ComplexF64}(undef,nkx,nky)
uekpar_star = Array{ComplexF64}(undef,nkx,nky)
gk_star = Array{ComplexF64}(undef,nkx,nky,ngtot)
dxg_star = Array{Float64}(undef,nlx,nly,ngtot)
dyg_star = Array{Float64}(undef,nlx,nly,ngtot)

akpar_new = Array{ComplexF64}(undef,nkx,nky)
uekpar_new = Array{ComplexF64}(undef,nkx,nky)
nek_new = Array{ComplexF64}(undef,nkx,nky)
phik_new = Array{ComplexF64}(undef,nkx,nky)
gk_new = Array{ComplexF64}(undef,nkx,nky,ngtot)
fgm_pred = Array{ComplexF64}(undef,nkx,nky,ngtot)

rel_error_array = Array{ComplexF64}(undef,nkx,nky)

if debugging
    print("Finished Initializing arrays \n")
end

# Begin setup of necessary values, dx,dys, hypercoeffs, timestep, Initial conditions

# Get initial conditions, load actual values into arrays!

apar_eq, phi_eq = equilibrium()

akpar_eq = FFT2d_direct(apar_eq) # will be used later!

#apar_perturb = init_perturb(apar_perturb) # TODO:Implement perturbation to phi or Apar
# apar = apar_eq + apar_perturb
# phi = phi_eq + phi_perturb
apar = apar_eq
phi = phi_eq 

akpar = FFT2d_direct(apar)
phik = FFT2d_direct(phi)

if rhoi < small_rhoi
    for i = 1:nkx
        for j = 1:nky
            nek_perturb[i,j] = -kperp(i,j)^2*phik[i,j]
        end
    end
else
    for i = 1:nkx
        for j = 1:nky
            nek_perturb[i,j] = 2.0/rhoi^2*(Γ₀(kperp(i,j)^2*rhoi^2/2.0)-1)*phik[i,j]
        end
    end
end
nek = nek_perturb

for i = 1:nkx
    for j = 1:nky
        uekpar[i,j] = -kperp(i,j)^2*akpar[i,j]
    end
end

# In Viriato, take invFFT to get ne, uepar in real space. But these only used for diagnostics, so dont do that here yet


# Take derivaties necesary to calculate hypercoeffs, timestep for first iteration

# TODO: MAKE SURE THE CONVOL DEFINITIONS ARE NOT CONVOL3,CONVOL4, just convol. Only need 2D! Im sure this will throw an error bc convol isnt defined!
dxphi,dyphi = convol(phik) # Better to write functions like this or as something that modifies value?
dxne, dyne = convol(nek)
dxapar,dyapar = convol(akpar)
dxuepar,dyuepar = convol(uekpar)

if g_inc 
    for ng = 1:ngtot
        dxg[:,:,ng],dyg[:,:,ng] = convol(gk[:,:,ng])
    end
end 

vex,vey,vrhosx,vrhosy = flows(dxphi,dyphi,dxne, dyne)
vxmax = max(maximum(abs.(vex)),maximum(abs.(vrhosx)))
vymax = max(maximum(abs.(vey)),maximum(abs.(vrhosy)))

bx,by = bfield(dxapar,dyapar)
bxmax = maximum(abs.(bx))
bymax = maximum(abs.(by))
bperp = sqrt.(bx.^2+by.^2)
bperp_max = maximum(bperp)

omega_kaw = omegakaw(bperp_max) # Function omega_kaw--its value is public in the function def

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
        η2=hyper_coef/dti/kperpmax^(2*hyper_order)
    end
end

if hyper_colls_fixed
    hyper_νei=hyper_colls
else    
    hyper_νei=hyperm_coef/dti/(ngtot+1)^(2*hyper_morder)
end

if debugging
    print("Calced hypercoeffs, ready to begin timeloop \n")
end

savetime = 0.0
############## TIME LOOP ###################

#First some values which must be initialized for LOOP
relative_error = 0.0
aa0 = init_aa0_fac 

#p = 0
#z = 0 # not necessary if 2d? Not sure if this is z as in z direction
repeat = false # Flags for certain behaviors in loop. Trying again with smaller timestep, for example
noinc = false
divergent = false
first = true # in priciple this should not be true if restarts are enabled, but lets get to that later
for t = 0:tmax
    #global repeat, noinc, divergent, first, phik, uekpar, nek, akpar,dti # to make sure julia doesn't make them local
    #p = t-z # ? Not sure that this is necessary, controls when some files are written for diagnostics...
    if repeat
        repeat = false
        t -= 1 # Want to redo the same timestep, so just reduce the t index by one. 
    else
        if divergent
            divergent = false
            t -= 1 # Want to redo the same timestep, so just reduce the t index by one.
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
                fgm_old[:,:,gmin] = func_g2(dxg[:,:,gmin],dyg[:,:,gmin],dxphi,dyphi,dxapar,dyapar,
                dxg[:,:,gmin+1],dyg[:,:,gmin+1],bracket_akpar_uekpar)
                # \mathcal{gm}
                for ng = gmin+1:ngtot-1
                    fgm_old[:,:,ng] .= func_gm(ng,dxg[:,:,ng-1],dyg[:,:,ng-1],dxg[:,:,ng],dyg[:,:,ng],dxg[:,:,ng+1],dyg[:,:,ng+1],
                    dxphi,dyphi,dxapar,dyapar) 
                end
                # \mathcal{glast}
                fglast_old = func_lastg(hyper_νei,η2,dxg[:,:,ngtot-1],dyg[:,:,ngtot-1],dxg[:,:,ngtot],dyg[:,:,ngtot],
                dxphi,dyphi,dxapar,dyapar)

            else 
                fApar_old = func_Akpar(dxapar,dyapar,dxphi,dyphi,dxne,dyne,dxuepar,dyuepar,dummy_real,dummy_real)
            end
        end 
    end

    # Get SI operator necessary for corrector step loop
    semi_implicit_operator = func_semi_implicit_operator(dti,bperp_max,aa0)

    
    # Start Predictor ("star") step

    if debugging
        print("Starting predictor step \n")
    end

    guess = akpar
    
    for i = 1:nkx
        for j = 1:nky
            nek_star[i,j] = exp_nu(i,j,ν2,dti)*nek[i,j]+ dti/2.0*(1.0+exp_nu(i,j,ν2,dti))*fne_old[i,j]

            akpar_star[i,j] = exp_eta(i,j,η2,dti)*akpar[i,j] + dti/2.0*(1.0+exp_eta(i,j,η2,dti))*fApar_old[i,j] + 
                (1.0-exp_eta(i,j,η2,dti))*akpar_eq[i,j]

            uekpar_star[i,j] = -kperp(i,j)^2*akpar_star[i,j]
        end
    end
    
    if g_inc
        # Get first and last g
        for i = 1:nkx
            for j = 1:nky
                gk_star[i,j,gmin] = exp_nu(i,j,ν2,dti)*gk[i,j,gmin] + dti/2.0*(1+exp_nu(i,j,ν2,dti))*fgm_old[i,j,gmin]
                
                gk_star[i,j,ngtot] = exp_ng(ngtot,hyper_νei,dti)*exp_nu(i,j,ν_g,dti)*gk[i,j,ngtot]+
                    dti/2.0*(1.0+exp_ng(ngtot,hyper_νei,dti)*exp_nu(i,j,ν_g,dti))*fglast_old[i,j]
            end 
        end
        # get the rest of the gs
        for ng = gmin+1:ngtot-1
            for i = 1:nkx
                for j = 1:nky
                    gk_star[i,j,ng] = exp_ng(ng,hyper_νei,dti)*exp_nu(i,j,ν_g,dti)*gk[i,j,ng]+
                        dti/2.0*(1.0+exp_ng(ng,hyper_νei,dti)*exp_nu(i,j,ν_g,dti))*fgm_old[i,j,ng]
                end
            end
        end

    end
    for ng = gmin:ngtot
        dxg_star[:,:,ng],dyg_star[:,:,ng] = convol(gk_star[:,:,ng])
    end 
    

    phik_star = phi_pot(nek_star)
    dxphi_star,dyphi_star = convol(phik_star) 
    dxne_star, dyne_star = convol(nek_star)
    dxapar_star,dyapar_star = convol(akpar_star)
    dxuepar_star,dyuepar_star = convol(uekpar_star)
    
    if g_inc
        fapar_star = func_Akpar(dxapar_star,dyapar_star,dxphi_star,dyphi_star,
            dxne_star,dyne_star,dxuepar_star,dyuepar_star,dxg_star[:,:,gmin],dyg_star[:,:,gmin])
    else
        fapar_star = func_Akpar(dxapar_star,dyapar_star,dxphi_star,dyphi_star,
            dxne_star,dyne_star,dxuepar_star,dyuepar_star,dummy_real,dummy_real)
    end
    
    

    # begin corrector "loop", although here we only do pmax = 1, so just do this once. Will lower timestep if not converged in timestep
    # if not converged after one step
    
    if debugging
        print("Starting corrector loop \n")
    end

    @debug "Starting corrector loop!!!!!!!!!"

    p_iter = 0
    for p_iter = 0:1
        p_count +=1
        sum_apar_rel_error = 0.0
        rel_error_array .= 0.0 # array of zeros?
        old_error = 0.0
        if p_iter == 0 # use star values as p=0
            fapar_pred = fapar_star
            value_nex = dxne_star
            value_ney = dyne_star
            value_phix = dxphi_star
            value_phiy = dyphi_star
            value_gx = dxg_star
            value_gy = dyg_star
        else
            old_error = relative_error # Old error stores relative error from last update. Use to check if divergent
            if g_inc
                fapar_pred = func_Akpar(dxapar,dyapar,dxphi,dyphi,dxne,dyne,dxuepar,dyuepar,dxg[:,:,gmin],dyg[:,:,gmin])
            else
                fapar_pred = func_Akpar(dxapar,dyapar,dxphi,dyphi,dxne,dyne,dxuepar,dyuepar,dummy_real,dummy_real)
            end

            # update values, RHS are calculated at bottom of this loop, ie from step p=0
            value_nex = dxne
            value_ney = dyne
            value_phix = dxphi
            value_phiy = dyphi
            value_gx = dxg
            value_gy = dyg
        end

        relative_error = 0.0 # reset value of relative error
        for i = 1:nkx
            for j = 1:nky
                akpar_new[i,j] = 1.0/(1.0+semi_implicit_operator[i,j]/4.0)*(exp_eta(i,j,η2,dti)*akpar[i,j]+
                    dti/2.0*exp_eta(i,j,η2,dti)*fApar_old[i,j]+ 
                    dti/2.0*fapar_pred[i,j]+
                    (1.0-exp_eta(i,j,η2,dti))*akpar_eq[i,j]+
                    semi_implicit_operator[i,j]/4.0*guess[i,j])
                
                uekpar_new[i,j] = -kperp(i,j)^2*akpar_new[i,j]
                
                # Take difference between akpar t=n and akpar_new, ie at t=n+1. Note akpar is only updated after p_loop. So this is akpar, n+1,p+1 - akpar,n
                sum_apar_rel_error = sum_apar_rel_error + (abs(akpar_new[i,j]-akpar[i,j]))^2 # Add up the change in Apar at every index. Will later divide by nkx*nky to get avg     
            end
        end

        dxapar,dyapar = convol(akpar_new)
        dxuepar,dyuepar = convol(uekpar_new)
        
        # to get next ne, use the last timetep values of everything except the new apar
        fne_pred, bracket_akpar_uekpar = func_ne(value_phix,value_phiy,value_nex,value_ney,dxapar,dyapar,dxuepar,dyuepar)
        for i = 1:nkx
            for j = 1:nky
                nek_new[i,j] = exp_nu(i,j,ν2,dti)*nek[i,j]+ dti/2.0*(1.0+exp_nu(i,j,ν2,dti))*fne_old[i,j] + dti/2.0*fne_pred[i,j]

                # Error for this p iteration at each location. Note that guess hold the value of akpar_new from previous p_loop iteration so this is akpar,n+1,p+1 - akpar,n+1,p 
                rel_error_array[i,j] = abs(semi_implicit_operator[i,j]/4.0*(akpar_new[i,j]-guess[i,j]))/sqrt(sum_apar_rel_error/(nkx*nky))
            end
        end
        
        # Get relative error for this p loop. Should definitely also put loop break statements here--no need to recalc g if we are only doing one p step!
        # If more than 1 p step, then would need to keep track of new g values over multiple p iterations...
        relative_error = maximum(abs.(rel_error_array))

        phik_new = phi_pot(nek_new)

        # update phi, ne derivatives
        dxphi, dyphi = convol(phik_new)
        dxne,dyne = convol(nek_new)

        
        # now get next mathcal gs 
        if g_inc
            
            #get g2 p+1
            fgm_pred[:,:,gmin] = func_g2(value_gx[:,:,gmin],value_gy[:,:,gmin],dxphi,dyphi,dxapar,dyapar,
            value_gx[:,:,gmin+1],value_gy[:,:,gmin+1],phik_new)
            for i = 1:nkx
                for j = 1:nky
                    gk_new[i,j,gmin] = exp_nu(i,j,ν2,dti)*gk[i,j,gmin] +
                        dti/2.0*exp_nu(i,j,ν2,dti)*fgm_old[i,j,gmin] + 
                        +dti/2.0*fgm_pred[i,j,gmin]
                end
            end

            dxg[:,:,gmin],dyg[:,:,gmin] = convol(gk_new[:,:,gmin])
            # get gm p+1
            for ng = gmin+1:ngtot-1
                
                dxgm = dxg[:,:,ng-1] # need this bc gm p+1 relies on gm-1 p+1
                dygm = dyg[:,:,ng-1]
                fgm_pred[:,:,ng] .= func_gm(ng,dxgm,dygm,value_gx[:,:,ng],value_gy[:,:,ng],value_gx[:,:,ng+1],value_gy[:,:,ng+1],
                dxphi,dyphi,dxapar,dyapar)

                for i = 1:nkx
                    for j = 1:nky 
                        gk_new[i,j,ng] = exp_ng(ng,hyper_νei,dti)*exp_nu(i,j,ν_g,dti)*gk[i,j,ng]+
                            dti/2.0*exp_ng(ng,hyper_νei,dti)*exp_nu(i,j,ν_g,dti)*fgm_old[i,j,ng]+
                            dti/2.0*fgm_pred[i,j,ng]
                    end
                end

                dxg[:,:,ng],dyg[:,:,ng] = convol(gk_new[:,:,ng])
            end
            
            # get glast p+1
            f_lastg = func_lastg(hyper_νei,η2,dxg[:,:,ngtot-1],dyg[:,:,ngtot-1],value_gx[:,:,ngtot],value_gy[:,:,ngtot],
                dxphi,dyphi,dxapar,dyapar)
            for i = 1:nkx
                for j = 1:nky
                    gk_new[i,j,ngtot] = exp_ng(ngtot,hyper_νei,dti)*exp_nu(i,j,ν_g,dti)*gk[i,j,ngtot]+
                        dti/2.0*exp_ng(ngtot,hyper_νei,dti)*exp_nu(i,j,ν_g,dti)*fglast_old[i,j]+
                            dti/2.0*f_lastg[i,j]
                end
            end

            dxg[:,:,ngtot],dyg[:,:,ngtot] = convol(gk_new[:,:,ngtot])
            
        end    
            
        # Now have all necessary values at p=1

        # Test for convergence
        if p_iter >= 1 && relative_error/old_error >= 1.0
            dti = low*dti
            #z=z+1 # maybe not actually necessary.
            divergent = true
            # exit ploop --> how to do in julia?
            break
        end

        if relative_error <= epsilon
            break # exit p loop
        end

        if relative_error > epsilon && p_iter==pmax
            dti=low*dti
            #z=z+1
            repeat= true
            break # exit p loop
        end

        guess = akpar_new
    end # end of p loop

    if debugging
        print("End of corrector loop \n")
    end

    if divergent
        continue # go to next time loop iteration with divergent = true
    end
    
    if repeat
        redo_timestep = true
        noinc = true
        continue # go to next time loop iteraton with repeat and noinc true
    end 

    # Update Variables to "new" values (i.e. p+1)

    nek = nek_new
    akpar = akpar_new
    phik = phik_new
    uekpar = uekpar_new

    savetime = savetime + dti # update the simulation timestep

    # Now re-evaluate the flows to evaluate CFL condition
    
    vex,vey,vrhosx,vrhosy = flows(dxphi,dyphi,dxne,dyne)
    vxmax = max(maximum(abs.(vex)),maximum(abs.(vrhosx)))
    vymax = max(maximum(abs.(vey)),maximum(abs.(vrhosy)))

    bx,by = bfield(dxapar,dyapar)
    bxmax = maximum(abs.(bx))
    bymax = maximum(abs.(by))
    bperp = sqrt.(bx.^2+by.^2)
    bperp_max = maximum(bperp)

    omega_kaw = omegakaw(bperp_max) # Function omega_kaw--its value is public in the function def

    # Seems to be just a diagnostic.
    # uxavg = (uxavg+ vex)/(p+1.0)
    # uyavg = (uyave + vey)/(p+1.0)

    # Calculate CFL fraction/timestep. Not sure why so simple here...
    CFL_flow=min(dx/vxmax,dy/vymax)

    dti_temp = CFL_frac*CFL_flow # TIMESTEP!

    # Here can call diagnostics if necessary, do I/O stuff
    
    # Calculate new timestep!
    dti = dtnext(relative_error,dti_temp,noinc,dti) 

    # Calc new hyper coeffs
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
            η2=hyper_coef/dti/kperpmax^(2*hyper_order)
        end
    end

    if hyper_colls_fixed
        hyper_νei=hyper_colls
    else    
        hyper_νei=hyperm_coef/dti/(ngtot+1)^(2*hyper_morder)
    end

    if debugging
        print("End of timeloop, moving to next iteration \n")
        print("Timestep ", t, "\n")
    end
end # END OF TIMELOOP

if debugging
    print("End of REGK \n")
end
close(log_file)
end # End of main 

main()