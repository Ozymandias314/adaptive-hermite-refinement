using LinearAlgebra, FFTW, Logging, JLD2

include("constants.jl")


function energy_tot(akpar,phik)
    b_energy_tot = 0.0
    phine_energy_tot = 0.0
    for i = 1:nkx
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
    return b_energy_tot,phine_energy_tot
end 

#phine_energy(k)-1.0/rhoi**2*(gama0(kperp(j,1)**2*rhoi**2/2.)-1.)*&
#FIK(j,1,k)*conjg(FIK(j,1,k))
