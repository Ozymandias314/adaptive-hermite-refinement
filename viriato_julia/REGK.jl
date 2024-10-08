using LinearAlgebra, FFTW, Logging, JLD2, Printf

function print_cpp(arr::Matrix{ComplexF64})
    for i in 1:size(arr, 1)
        row_str = ""
        for j in 1:size(arr, 2)
            real_part = real(arr[i, j])
            imag_part = imag(arr[i, j])
            row_str *= "($(real_part), $(imag_part)) "
        end
        println(row_str)
    end
end

function print_cpp(arr::Matrix{Float64})
    for i in 1:size(arr, 1)
        row_str = ""
        for j in 1:size(arr, 2)
            row_str *= "$(arr[i, j]) "
        end
        println(row_str)
    end
end


include("constants.jl")
include("transforms.jl")
include("Brackets.jl")
include("aux.jl")
include("diag.jl")
include("functions.jl")
include("initialize.jl")
include("diagnostics.jl") # Right now this just calculates energy
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


debugging = false # Debugging print statements


phik = Array{ComplexF64}(undef,nkx,nky)
nek = Array{ComplexF64}(undef,nkx,nky)
akpar = Array{ComplexF64}(undef,nkx,nky)
uekpar = Array{ComplexF64}(undef,nkx,nky)
nek_perturb = Array{ComplexF64}(undef,nkx,nky)

dummy_real = zeros(Float64,nlx,nly)

gk = zeros(ComplexF64,nkx,nky,ngtot)
dxg = Array{Float64}(undef,nlx,nly,ngtot)
dyg = Array{Float64}(undef,nlx,nly,ngtot)
fgm_old = Array{ComplexF64}(undef,nkx,nky,ngtot)
if debugging
    print("Gk","\n")
    print(gk[32,32,gmin],'\n')
end
nek_star = Array{ComplexF64}(undef,nkx,nky)
phik_star = Array{ComplexF64}(undef,nkx,nky)
akpar_star = Array{ComplexF64}(undef,nkx,nky)
uekpar_star = Array{ComplexF64}(undef,nkx,nky)
gk_star = zeros(ComplexF64,nkx,nky,ngtot)
dxg_star = Array{Float64}(undef,nlx,nly,ngtot)
dyg_star = Array{Float64}(undef,nlx,nly,ngtot)

akpar_new = Array{ComplexF64}(undef,nkx,nky)
uekpar_new = Array{ComplexF64}(undef,nkx,nky)
nek_new = Array{ComplexF64}(undef,nkx,nky)
phik_new = Array{ComplexF64}(undef,nkx,nky)
gk_new = zeros(ComplexF64,nkx,nky,ngtot)
fgm_pred = Array{ComplexF64}(undef,nkx,nky,ngtot)

# Initialze these variables here, otherwise julia forgets about them between for loop iterations when Repeat or Divergent are activated
fne_old = undef
fApar_old = undef
fglast_old = undef


rel_error_array = Array{ComplexF64}(undef,nkx,nky)

if debugging
    print("Finished Initializing arrays \n")
end

# Begin setup of necessary values, dx,dys, hypercoeffs, timestep, Initial conditions

# Get initial conditions, load actual values into arrays!

apar_eq, phi_eq = equilibrium()

akpar_eq = FFT2d_direct(apar_eq) # will be used later!

apar_perturb = init_perturb() # TODO:Implement perturbation to phi or Apar
apar = apar_eq + apar_perturb
phi = phi_eq

file_string_apar = "apar_initial.jld2"
file_string_ne = "ne_initial.jld2"
# Save Apar, ne in real space
#save_object(file_string_apar,apar)


if debugging
    print("initialized values","\n")
    print(apar_eq[32,32],apar[32,32],phi[32,32],'\n')
    println("Initial maximum akpar", maximum(abs.(akpar)))
    println()
    println("Initial maximum apar", maximum(abs.(apar)))
    println()
    println("Initial maximum apar from iFFT akpar ", maximum(abs.(FFT2d_inv(akpar))))
end

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

if debugging
    print("Data from initial convol call","\n")
    print(dxuepar[32,32],dxapar[32,32],dxne[32,32],dxuepar[32,32],'\n')
    print(gk[32,32,gmin],"\n")
end

if g_inc 
    for ng = 1:ngtot
        dxg[:,:,ng],dyg[:,:,ng] = convol(gk[:,:,ng])
    end
end 

if debugging
    print("initial dxg dyg calls","\n")
    print(dxg[32,32,gmin],dyg[32,32,gmin],gk[32,32,gmin],'\n')
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
# TODO Viriato has terms with dz here, which I'm not sure are actually zero. Make sure!
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
low = 0.92
############## TIME LOOP ###################

#First some values which must be initialized for LOOP
relative_error = 0.0
aa0 = init_aa0_fac 

#p = 0
#z = 0 # not necessary if 2d? Not sure if this is z as in z direction
repeat = false # Flags for certain behaviors in loop. Trying again with smaller timestep, for example
count_repeats = 0
count_divergent = 0
noinc = false
divergent = false
first = true # in priciple this should not be true if restarts are enabled, but lets get to that later
t = 0
while t <= tmax
    #global repeat, noinc, divergent, first, phik, uekpar, nek, akpar,dti # to make sure julia doesn't make them local
    #p = t-z # ? Not sure that this is necessary, controls when some files are written for diagnostics...
    if repeat
        repeat = false
        count_repeats += 1
        #t -= 1 # Want to redo the same timestep, so just reduce the t index by one. 
    else
        if divergent
            #print("Here divergent true")
            divergent = false
            count_divergent += 1
            #t -= 1 # Want to redo the same timestep, so just reduce the t index by one.
        else

            #p_count = 0 # Number of loops through corrector step 
            
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

                first = false
            end
            

            # Get the nonlinear operator values
            
            #\mathcal{N}
            fne_old, bracket_akpar_uekpar = func_ne(dxphi,dyphi,dxne,dyne,dxapar,dyapar,dxuepar,dyuepar)
            # "old" values are calculated from the previous timestep. These values, dxphi,dyphi, etc are only updated after succesful timestep

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

    if debugging
        print("Data from first func cals","\n")
        print(fApar_old[32,32],fgm_old[32,32,1],fne_old[32,32],'\n')

        print("Inputs to fApar","\n")
        print(dxapar[32,32],dyapar[32,32,1],dxphi[32,32],dyphi[32,32],dxne[32,32],dyne[32,32],dxuepar[32,32],dyuepar[32,32],dxg[32,32,gmin],dyg[32,32,gmin],'\n')
    end

    # Get SI operator necessary for corrector step loop
    semi_implicit_operator = func_semi_implicit_operator(dti,bperp_max,aa0)
    # Start Predictor ("star") step. Star values are the update from the predictor step, calculated based on the previous timestep. 

    if debugging
        print("Starting predictor step \n")
    end

    guess = deepcopy(akpar)

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
                gk_star[i,j,gmin] = exp_nu(i,j,ν2,dti)*gk[i,j,gmin] + dti/2.0*(1.0+exp_nu(i,j,ν2,dti))*fgm_old[i,j,gmin]
                
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
    
    if debugging
        print("Data from star step","\n")
        print(phik_star[32,32],nek_star[32,32],akpar_star[32,32],'\n')
    end

    # begin corrector "loop", although here we only do pmax = 1, so just do this once. Will lower timestep if not converged in timestep
    # if not converged after one step
    
    if debugging
        print("Starting corrector loop \n")
    end

    @debug "Starting corrector loop!!!!!!!!!"

    p_iter = 0
    for p_iter = 0:pmax
        #p_count +=1
        sum_apar_rel_error = 0.0
        rel_error_array .= 0.0 # array of zeros?
        old_error = 0.0
        if p_iter == 0 # use star values as p=0
            fapar_pred = fapar_star
            value_nex = dxne_star
            value_ney = dyne_star
            value_phix = dxphi_star
            value_phiy = dyphi_star
            if g_inc
                value_gx = dxg_star
                value_gy = dyg_star
            end
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
            if g_inc
                value_gx = dxg
                value_gy = dyg
            end
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
                sum_apar_rel_error = sum_apar_rel_error + abs2(akpar_new[i,j]-akpar[i,j]) # Add up the change in Apar at every index. Will later divide by nkx*nky to get avg     
            end
        end

        dxapar,dyapar = convol(akpar_new)
        dxuepar,dyuepar = convol(uekpar_new)
        
        if debugging
        println("maximum diff is ",maximum(abs.(akpar_new-akpar)))
        end

        # to get next ne, use the last timetep values of everything except the new apar
        fne_pred, bracket_akpar_uekpar = func_ne(value_phix,value_phiy,value_nex,value_ney,dxapar,dyapar,dxuepar,dyuepar)
        for i = 1:nkx
            for j = 1:nky
                nek_new[i,j] = exp_nu(i,j,ν2,dti)*nek[i,j]+ dti/2.0*exp_nu(i,j,ν2,dti)*fne_old[i,j] + dti/2.0*fne_pred[i,j]

                # Error for this p iteration at each location. Note that guess hold the value of akpar_new from previous p_loop iteration so this is akpar,n+1,p+1 - akpar,n+1,p 
                rel_error_array[i,j] = abs(semi_implicit_operator[i,j]/4.0*(akpar_new[i,j]-guess[i,j]))/sqrt(sum_apar_rel_error/(nkx*nky))
            end
        end
        
        # Get relative error for this p loop. Should definitely also put loop break statements here--no need to recalc g if we are only doing one p step!
        # If more than 1 p step, then would need to keep track of new g values over multiple p iterations...
        relative_error = maximum(abs.(rel_error_array))

#         if debugging

        println("sumApar relative_error: $sum_apar_rel_error")
        println("relative_error: $relative_error")
#         end

        phik_new = phi_pot(nek_new)

        # update phi, ne derivatives
        dxphi, dyphi = convol(phik_new)
        dxne,dyne = convol(nek_new)

        
        # now get next mathcal gs 
        if g_inc
            
            #get g2 p+1
            fgm_pred[:,:,gmin] = func_g2(value_gx[:,:,gmin],value_gy[:,:,gmin],dxphi,dyphi,dxapar,dyapar,
            value_gx[:,:,gmin+1],value_gy[:,:,gmin+1],bracket_akpar_uekpar)
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
        guess = deepcopy(akpar_new) # guess only updated if error is low enough from this iteration. 
    end # end of p loop

    if debugging
        print("End of corrector loop \n")
    end

    if divergent
        continue # go to next time loop iteration with divergent = true
    end
    
    if repeat
        noinc = true # Tells the timestep update if the timestep is allowed to increase
        continue # go to next time loop iteraton with repeat and noinc true
    end 

    # Update Variables to "new" values (i.e. p+1)

    if debugging
        print("Data at new timestep","\n")
        print(phik_new[32,32],nek_new[32,32],akpar_new[32,32],'\n')
    end

    nek = deepcopy(nek_new)
    akpar = deepcopy(akpar_new)
    phik = deepcopy(phik_new)
    uekpar = deepcopy(uekpar_new)
    if g_inc 
        gk = deepcopy(gk_new)
    end
    savetime += dti # update the simulation timestep

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

    # Calculate CFL fraction/timestep. Not sure why so simple here, but this is how its done for the not "turb" case in Viriato
    CFL_flow=min(dx/vxmax,dy/vymax)

    dti_temp = CFL_frac*CFL_flow # TIMESTEP!
    
    # Here can call diagnostics if necessary, do I/O stuff
    
    # Calculate new timestep!
    dti,noinc = dtnext(relative_error,dti_temp,noinc,dti) 
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

    b_energy_tot,phine_energy_tot = energy_tot(akpar,phik)
    println("magnetic energy: $b_energy_tot, kinetic energy: $phine_energy_tot")

    println("----------")
    println("Moving on to next timestep, ", t+1)
    println("dti is ",dti)
    #println("savetime is ", savetime)
    println("num repeats is ", count_repeats)
    println("Divergent counts ",count_divergent)
    println("----------")
    count_repeats=0
    count_divergent=0

    if debugging
        print("Final Data","\n")
        print(phi[32,32],nek[32,32],akpar[32,32],'\n')
    end
    #print("At end of tloop, t = ",t)
    # DIAGNOSTICS GO HERE
    if save_datafiles != 0 && t%save_datafiles == 0
        file_string_apar = "apar_"*string(t)*".jld2"
        #file_string_ne = "ne_"*string(t)*".jld2"
        # Save Apar, ne in real space
        apar = FFT2d_inv(akpar)
        #ne = FFT2d_inv(nek)
        save_object(file_string_apar,apar)
        #save_object(file_string_ne,ne)
        println("Saved data for timestep = ",t, " savetime= ",savetime)
        # println("relative_error ", relative_error) 
        # println("dti ",dti," temp dti ", dti_temp) # Factor of 2 small for dti_temp-->direct calc from flows . Factor of 5.7 small for actual timesteps --> why? Has to do w relative error as well, but relative error very similar
        # println("bxmax,bymax,bperpmax ",bxmax," ",bymax," ",bperp_max )
    end


    t += 1
end # END OF TIMELOOP


if debugging
    print("End of REGK \n")
    println("Repeat counts ",count_repeats)
    println("Divergent counts ",count_divergent)
end
close(log_file)
end # End of main 

main()