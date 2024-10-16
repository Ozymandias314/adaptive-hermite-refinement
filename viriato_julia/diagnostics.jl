using LinearAlgebra, FFTW, Logging, JLD2

include("constants.jl")


function energy_tot(akpar,phik)
    b_energy_tot = 0.0
    phine_energy_tot = 0.0
    
    for j=1:nky
        b_energy_tot += 0.5*kperp(1,j)^2*akpar[1,j]*conj(akpar[1,j])
        if rhoi < small_rhoi
            phine_energy_tot += 1.0*kperp(1,j)^2*phik[1,j]*conj(phik[1,j])
        else
            phine_energy_tot -= 1.0/rhoi^2*(Γ₀(kperp(1,j)^2*rhoi^2/2.0)-1)*phik[1,j]*conj(phik[1,j])
        end
    end
    
    for i = 2:nkx
        for j = 1:nky
            # Technically not correct, need to deal with double counting symmetric mode
            # Julia retains some complex part after this multiplication...numerical issue?
            b_energy_tot += kperp(i,j)^2 * abs2(akpar[i,j])
            if rhoi < small_rhoi
                phine_energy_tot += 1.0*kperp(i,j)^2 * abs2(phik[i,j])
            else
                phine_energy_tot -= 1.0/rhoi^2*(Γ₀(kperp(i,j)^2*rhoi^2/2.0)-1) * abs2(phik[i,j])
            end
        end 
    end 
    return b_energy_tot/nlx/nly,phine_energy_tot/nlx/nly
end 

#phine_energy(k)-1.0/rhoi**2*(gama0(kperp(j,1)**2*rhoi**2/2.)-1.)*&
#FIK(j,1,k)*conjg(FIK(j,1,k))

# For now just the magnetic, electric, kinetic energies
function k_energy(akpar,phik)
    
    energy_u = Array{Float64}(undef,nkx,nky)
    energy_e = Array{Float64}(undef,nkx,nky)
    energy_b = Array{Float64}(undef,nkx,nky)
    kshell_u = zeros(Float64,ceil(kperpmax))
    kshell_e = zeros(Float64,ceil(kperpmax))
    kshell_b = zeros(Float64,ceil(kperpmax))
    kp_array = collect(range(1,ceil(kperpmax)))

    
    # First calculate the energies in k-space

    for j = 1:nky
        energy_b[1,j] = 0.5*abs2(kperp(1,j)*akpar[1,j])
        energy_e[1,j] = 0.5*abs2(kperp(1,j)*phik[1,j])
        energy_u[1,j] = -1.0/rhoi^2*(Γ₀(kperp(1,j)^2*rhoi^2/2.0)-1) * abs2(phik[1,j])
    end

    for i = 2:nkx
        for j = 1:nky

            energy_b[i,j] = abs2(kperp(i,j)*akpar[i,j])
            energy_e[i,j] = abs2(kperp(1,j)*phik[1,j])
            energy_u[i,j] = -2.0/rhoi^2*(Γ₀(kperp(1,j)^2*rhoi^2/2.0)-1) * abs2(phik[1,j])

            
        end 
    end 

    # Now get the perpendicular spectrum

    # Viriato version 
    # for kp = 1:ceil(kperpmax)
    #     for i = 1:nkx
    #         for j = 1:nky

    #             if (sqrt((i-1)^2+ky(j)^2) >= kp-1 && sqrt((i-1)^2+ky(j)^2) < kp+1 )
    #                 kshell_b(kp) += energy_b[i,j]
    #                 kshell_e(kp) += energy_e[i,j]
    #                 kshell_u(kp) += energy_u[i,j]
    #             end
    #         end
    #     end
    # end

    # Faster version
    for i = 1:nkx
        for j = 1:nky
            kp = round(Int,(sqrt((i-1)^2+ky(j)^2)))
            kshell_b[kp] += energy_b[i,j]
            kshell_e[kp] += energy_e[i,j]
            kshell_u[kp] += energy_u[i,j]
        end
    end
    
    return kp_array,kshell_b,kshell_e,kshell_u # A list of the kp values and the energies contained in each shell
end

function gm_spectrum(gk)

    energy_gm = zeros(Float64,(ngtot-gmin))
    gm_array = collect(range(1,(ngtot-gmin)))

    for m = gmin:ngtot

        for j = 1:nky
            energy_gm[m] += 0.5*abs2(rhos_diag*gk[1,j,m])
        end
        for i = 2:nkx
            for j = 1:nky
                energy_gm[m] += abs2(rhos_diag*gk[i,j,m])
            end
        end 

    end
    
    return gm_array, energy_gm # A list of the moment numbers and the energy contained in each
end