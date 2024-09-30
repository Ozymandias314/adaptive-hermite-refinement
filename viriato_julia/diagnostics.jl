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
