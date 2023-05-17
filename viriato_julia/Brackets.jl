#  use Brackets,   only: funcne_i, funcAkpar_i, funcg2, funcgm, func_lastg, bracket_3
#print(lambda)

# Only need bracket 3. Also, dont need the nlz stuff, as we are only
function Bracket_3(dxF, dyF, dxG, dyG)
    # this function calculates the "bracket" of two quantities, i.e., 
    # dxk/dx*dyk/dy-dxk/dy*dyk/dx. inputs are xk and yk, the quantities in
    # k space. Returns braxyk, the value of the bracket also in k-space.

    braxyk = zeros(nkx, nky)
    braxy = Array{Float64}(undef, nlx, nly)

    @. braxy = dxF * dyG - dyF * dxG

    braxyk = FFT2d_direct(braxy)

    braxyk[1, 1] = 0.0 # ensures that no zeroth mode is created

    return braxyk
end

#ok we need to figure out how the Bracket interface in FORTRAN is overloaded. I still don't get it. 

# Calculates nonlinear operator in 0th moment equation (Eq 44)
function func_ne(Dxphi, Dyphi, Dxne, Dyne, DxApar, DyApar, Dxuepar, Dyuepar)
    fne = zeros(ComplexF64, nkx, nky)

    #ya need to figure out bracket!
    bracket_phik_nek = Bracket_3(Dxphi, Dyphi, Dxne, Dyne)
    bracket_akpar_uekpar = Bracket_3(DxApar, DyApar, Dxuepar, Dyuepar)

    @. fne = -bracket_phik_nek + bracket_akpar_uekpar

    return fne, bracket_akpar_uekpar
end

# Calculates nonlinear operator in 1st moment equation (Eq 45)
function func_Akpar(DxApar,DyApar,Dxphi,Dyphi,Dxne,Dyne,Dxuepar,Dyuepar,Dxg2,Dyg2)
    fapar = zeros(ComplexF64, nkx, nky)

    tempx = Dxphi - rhos^2*(Dxne +sqrt(2.0)*Dxg2 )
    tempy = Dyphi - rhos^2*(Dyne +sqrt(2.0)*Dyg2 )

    bracket_akpar_phik = Bracket_3(DxApar, DyApar, tempx, tempy)
    bracket_uekpar_phik = Bracket_3(Dxuepar, Dyuepar, Dxphi, Dyphi)

    # these lines have like some random constants that are never defined? or commente dout? idk whats going on there 
    # also braakparuekpar never used in this subroutine (yes should be brauekparphi FIXED). For our purposes, notanj = 1.0 and rhoe_LTe = 0.0
    for j = 1:nky
        for i = 1:nkx
            # fapar[j, i, k] = 1.0/(1.0+kperp[j, i]^2*de^2)*
            #     (braakparphik[j, i, k] - de^2*brauekparphik[j, i, k] -
            #      notanj*1.0/sqrt(2.0)*rhos_de*rhoe_LTe*(0.0,1.0)*ky[j]*akpar[j, i, k])
            fapar[i, j] = 1.0/(1.0+kperp(i,j)^2*de^2)*
                (bracket_akpar_phik[i, j] - de^2*bracket_uekpar_phik[i, j])
        end
    end


    return fapar

end 
# Calculates nonlinear operator in 2nd moment
function func_g2(Dxg2, Dyg2, Dxphi, Dyphi, Dxapar, Dyapar, Dxg3, Dyg3, bracket_akpar_uekpar)
    nlx, nly = size(Dxapar)
    nkx, nky = size(bracket_akpar_uekpar)

    fg2 = zeros(ComplexF64, nkx, nky)

    bracket_g2_phik = Bracket_3(Dxg2, Dyg2, Dxphi, Dyphi)
    bracket_akpar_g3 = Bracket_3(Dxapar, Dyapar, Dxg3, Dyg3)

    for j in 1:nky, i in 1:nkx
        fg2[i,j] = bracket_g2_phik[i,j] + sqrt(gmin+1.0) * rhos_de * bracket_akpar_g3[i,j] +
                     notanj * sqrt(2.0) * bracket_akpar_uekpar[i,j] #-
                     #notanj * 1.0/(2.0*rhos_de) * rhoe_LTe*(0.0,1.0) * ky[j] * phik[j,i,k]
    end

    return fg2
end

#commented out the below block because it seems like it is using modules, and maybe we can just use include? or if we want to use module
# have to actually export correctly
#=
using ..constants: gmin, ngtot, npe, nlx, nly, nly_par, nlz_par, nky, nkx_par,
                   kpar0, j1, omega0, facpm, amplitude, Lz, pi, lambda, de, notanj,
                   rhos_de, rhoe_lte
using ..grid: ky, zz
using ..mp: iproc, proc0
using ..Functions: anj_kron
=#
function func_gm(m, Dxgm, Dygm, Dxg, Dyg, Dxgp, Dygp, Dxphi, Dyphi, Dxapar, Dyapar)

    fgm = zeros(ComplexF64, nkx, nky)
    #... Local vars
    bracket_akpar_gpm = zeros(ComplexF64, nkx, nky)
    bracket_phik_g = zeros(ComplexF64, nkx, nky)
    

    #feel like this part might be wrong....
    lte_kron = zeros(Float64, ngtot)
    lte_kron[gmin+1] = 1.0

    bracket_phik_g = Bracket_3(Dxphi, Dyphi, Dxg, Dyg)

    bracket_akpar_gpm = Bracket_3(Dxapar, Dyapar, 
        sqrt((m+1)*1.0) .* Dxgp .+ sqrt(m*1.0) .*  Dxgm, 
        sqrt((m+1)*1.0) .* Dygp .+ sqrt(m*1.0) .*  Dygm)

    for j = 1:nky
        for i = 1:nkx
            fgm[i, j] = -bracket_phik_g[i, j] + rhos_de * bracket_akpar_gpm[i, j]
        end
    end

    return fgm
end

# Calculates nonlinear part of last moment closure. Not explicitly in paper, but described in Eq 22
function func_lastg(hyper_nuei, niu2, Dxgm, Dygm, Dxg, Dyg, Dxphi, Dyphi, Dxapar, Dyapar)
    # Calculates g_M using the analytical nonlinear closure
    # As is, only does the perp part
    
    # Use statements
    # Not sure about nlx, nly, nkx, nlz_par, nky, ngot, nu_ei, kperp, niu, and hyper_morder// should be from constants
    # Define constants and modules if necessary (the grid stuff idk)
    # Omitted for brevity
    
    # Input:
    # hyper_nuei: real
    # niu2: real
    # Dxgm: real[nlx, nly, nlz_par]
    # Dygm: real[nlx, nly, nlz_par]
    # Dxg: real[nlx, nly, nlz_par]
    # Dyg: real[nlx, nly, nlz_par]
    # Dxphi: real[nlx, nly, nlz_par]
    # Dyphi: real[nlx, nly, nlz_par]
    # Dxapar: real[nlx, nly, nlz_par]
    # Dyapar: real[nlx, nly, nlz_par]
    
    # Output:
    # fglast: complex[nky, nkx, nlz_par]
    
    fglast = zeros(ComplexF64,nkx,nky)

    bracket_phik_gk = Bracket_3(Dxphi, Dyphi, Dxg, Dyg)
    bracket_akpar_gk = Bracket_3(Dxapar, Dyapar, Dxg, Dyg)

    
    for j = 1:nky
        for i = 1:nkx
            bracket_akpar_gk[i,j] = (rhos_de^2 * (ngtot+1)) / ((ngtot+1) * nu_ei +
                (ngtot+1)^(2*hyper_morder) * hyper_nuei +
                niu * kperp(i,j)^2 + niu2*kperp(i,j)^(2*hyper_order)) * bracket_akpar_gk[i,j]
        end
    end
    
    dx_bracket_akpar_gk,dy_bracket_akpar_gk = convol(bracket_akpar_gk)
    total_bracket = Bracket_3(Dxapar,Dyapar,dx_bracket_akpar_gk + rhos_de*sqrt(ngtot*1.0)*Dxgm,
        dy_bracket_akpar_gk + rhos_de*sqrt(ngtot*1.0)*Dygm)
    
    fglast = -bracket_phik_gk + total_bracket

    return fglast

end

