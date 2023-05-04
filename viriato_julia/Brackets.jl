#  use Brackets,   only: funcne_i, funcAkpar_i, funcg2, funcgm, func_lastg, bracket_3
include("constants.jl")
print(lambda)

# NOT USED?
# function Funcne_i(dxphi, dyphi, dxne, dyne, dxapar, dyapar, dxuepar, dyuepar, fne, braakparuekpar)
#     # assuming the necessary constants are defined in a module called "constants"
#     nkx = constants.nkx
#     nly = constants.nly
#     nlz_par = constants.nlz_par
#     nlx = constants.nlx
#     nly = constants.nly
#     nky = constants.nky

#     # I don't think you need to declare them like that in julia ^, once you do include(constants.jl) you can just use the name directly

#     braphiknek = Array{ComplexF64}(undef, nky, nkx, nlz_par)
#     for k in 1:nlz_par, i in 1:nkx, j in 1:nky
#         braphiknek[j, i, k] = dxphi[i, j, k] + im * dyphi[i, j, k] + dxne[i, j, k] + im * dyne[i, j, k]
#     end

#     braakparuekpar = Array{ComplexF64}(undef, nky, nkx, nlz_par)
#     for k in 1:nlz_par, i in 1:nkx, j in 1:nky
#         braakparuekpar[j, i, k] = dxapar[i, j, k] + im * dyapar[i, j, k] + dxuepar[i, j, k] + im * dyuepar[i, j, k]
#     end

#     for k in 1:nlz_par, i in 1:nkx, j in 1:nky
#         fne[j, i, k] = -braphiknek[j, i, k] + braakparuekpar[j, i, k] # assuming the commented line is not required
#     end
#     return nothing  # assuming the "braakparuekpar" output argument was modified in-place
# end

# Not needed
function Bracket_4(dxF, dyF, dxG, dyG)
    # this subroutine calculates the "bracket" of two quantities, i.e., 
    # dxk/dx*dyk/dy-dxk/dy*dyk/dx. inputs are xk and yk, the quantities in
    # k space. returns
    # braxyk, the value of the bracket also in k-space.

    #i think we can infer that it has the same shape
    braxy = similar(dxG)

    ngmin = 1

    #or ngtot I guess
    ngmax = size(dxG, 4)

    if linear
        braxy .= 0.0
    else 
        @. braxy = dxF * dyG - dyF * dxG

        #call the fft2d_direct on braxy, braxyk

        if mod(iproc, npe) == 0
            braxyk[1, 1, :, ngmin:ngmax] .= 0.0 # insures that no zeroth mode is created
        end
    end
    return braxyk
end

# Only need bracket 3. Also, dont need the nlz stuff, as we are only
function Bracket_3(dxF, dyF, dxG, dyG)
    # this function calculates the "bracket" of two quantities, i.e., 
    # dxk/dx*dyk/dy-dxk/dy*dyk/dx. inputs are xk and yk, the quantities in
    # k space. Returns braxyk, the value of the bracket also in k-space.

    braxyk = zeros(nky, nkx, nlz_par)

    if linear 
        return braxyk
    else
        @. braxy = dxF * dyG - dyF * dxG

        #some weird ifdef gasca2d, gasca3d thing, seems like it just affects how the FFT is calculated

        if mod(iproc, npe) == 0
            braxyk[1, 1, :, :] .= 0.0 # insures that no zeroth mode is created
        end
        return braxyk
    end
end

#ok we need to figure out how the Bracket interface in FORTRAN is overloaded. I still don't get it. 

# Calculates nonlinear operator in 0th moment equation (Eq 44)
function func_ne(Dxphi, Dyphi, Dxne, Dyne, DxApar, DyApar, Dxuepar, Dyuepar)
    fne = zeros(ComplexF64, nkx, nky)

    #ya need to figure out bracket!
    bracket_phik_nek = bracket(Dxphi, Dyphi, Dxne, Dyne)
    bracket_akpar_uekpar = bracket(DxApar, DyApar, Dxuepar, Dyuepar)

    @. fne = -bracket_phik_nek + bracket_akpar_uekpar

    return fne, bracket_akpar_uekpar
end

# Calculates nonlinear operator in 1st moment equation (Eq 45)
function func_Akpar(DxApar,DyApar,Dxphi,Dyphi,Dxne,Dyne,Dxuepar,Dyuepar,Dxg2,Dyg2)
    FApar = zeros(ComplexF64, nkx, nky)

    tempx = Dxphi - rhos^2*(Dxne +sqrt(2.0)*Dxg2 )
    tempy = Dyphi - rhos^2*(Dyne +sqrt(2.0)*Dyg2 )

    bracket_akpar_phik = bracket(DxApar, DyApar, tempx, tempy)
    bracket_uekpar_phik = bracket(Dxuepar, Dyuepar, Dxphi, Dyphi)

    # these lines have like some random constants that are never defined? or commente dout? idk whats going on there 
    # also braakparuekpar never used in this subroutine (yes should be brauekparphi FIXED). For our purposes, notanj = 1.0 and rhoe_LTe = 0.0
    for i = 1:nkx
        for j = 1:nky
            # fapar[j, i, k] = 1.0/(1.0+kperp[j, i]^2*de^2)*
            #     (braakparphik[j, i, k] - de^2*brauekparphik[j, i, k] -
            #      notanj*1.0/sqrt(2.0)*rhos_de*rhoe_LTe*(0.0,1.0)*ky[j]*akpar[j, i, k])
            fapar[i, j] = 1.0/(1.0+kperp[i, j]^2*de^2)*
                (bracket_akpar_phik[i, j] - de^2*bracket_uekpar_phik[i, j])
        end
    end


    return FApar

end 
# Calculates nonlinear operator in 2nd moment
function func_g2(Dxg2, Dyg2, Dxphi, Dyphi, Dxapar, Dyapar, Dxg3, Dyg3, bracket_akpar_uekpar)
    nlx, nly = size(Dxapar)
    nkx, nky = size(braakparuekpar)

    fg2 = zeros(ComplexF64, nkx, nky)

    bracket_g2_phik = bracket(Dxg2, Dyg2, Dxphi, Dyphi)
    bracket_akpar_g3 = bracket(Dxapar, Dyapar, Dxg3, Dyg3)

    for i in 1:nkx, j in 1:nky
        fg2[i,j] = bracket_g2_phik[i,j] + sqrt(gmin+1.0) * rhos_de * bracket_akpar_g3[i,j] +
                     notanj * sqrt(2.0) * bracket_akpar_uekpar[i,j] #-
                     #notanj * 1.0/(2.0*rhos_de) * rhoe_LTe*(0.0,1.0) * ky[j] * phik[j,i,k]
    end

    return fg2
end


using ..constants: gmin, ngtot, npe, nlx, nly, nly_par, nlz_par, nky, nkx_par,
                   kpar0, j1, omega0, facpm, amplitude, Lz, pi, lambda, de, notanj,
                   rhos_de, rhoe_lte
using ..grid: ky, zz
using ..mp: iproc, proc0
using ..Functions: anj_kron

function Funcgm_3(m, Dxgm, Dygm, Dxg, Dyg, Dxgp, Dygp, Dxphi, Dyphi, Dxapar, Dyapar, akpar, t)



    #...  intent OUT
    #...  size is nky * nkx_par * nlz_par
    fgm = zeros(ComplexF64, nkx_par, nky, nlz_par)
    #... Local vars
    braakpargpm = zeros(ComplexF64, nky, nkx_par, nlz_par)
    braphikg = zeros(ComplexF64, nky, nkx_par, nlz_par)
    Akpar = zeros(ComplexF64, nky, nkx_par, nlz_par)

    #feel like this part might be wrong....
    lte_kron = zeros(Float64, ngtot)
    lte_kron[gmin+1] = 1.0

    Bracket(Dxphi, Dyphi, Dxg, Dyg, braphikg)

    Bracket(DxApar, DyApar, sqrt((m+1)*1.0) .* Dxgp .+ sqrt(m*1.0) .* (1.0-1.0/lambda*anj_kron(m)*(1.0-notanj)) .* Dxgm, sqrt((m+1)*1.0) .* Dygp .+ sqrt(m*1.0) .* (1.0-1.0/lambda*anj_kron(m)*(1.0-notanj)) .* Dygm
    , braakpargpm)

    for k = 1:nlz_par
        for i = 1:nkx_par
            for j = 1:nky
                fgm[j, i, k] = -braphikg[j, i, k] + rhos_de * braakpargpm[j, i, k] +
                    notanj * lte_kron[m] * sqrt(3.0)/(2.0*de^2) *
                    rhoe_lte*(0.0, 1.0) * ky[j] * Akpar[j, i, k]
                # +notanj*lte_kron(m)*sqrt(3.0/2.0)*rhos_de*rhoe_LTe*(0.0,1.0)*ky(j)*Akpar(j,i,k)
                # NFL: 24/05/13 commented line above had wrong normalizations
            end
        end
    end
end

# ok now funcgm is used by two interfaces. not gonna do it now since I think they r related to the dimensionality and maybe we only need to do one. 

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

    bracket_phik_gk = Bracket(Dxphi, Dyphi, Dxg, Dyg)
    bracket_akpar_gk = Bracket(Dxapar, Dyapar, Dxg, Dyg)

    
    for i = 1:nkx
        for j = 1:nky
            bracket_akpar_gk[i,j] = (rhos_de^2 * (ngtot+1)) / ((ngtot+1) * nu_ei +
                (ngtot+1)^(2*hyper_morder) * hyper_nuei +
                niu * kperp[i,j]^2 + niu2*kperp[i,j]^(2*hyper_order)) * bracket_akpar_gk[i,j]
        end
    end
    
    dx_bracket_akpar_gk,dy_bracket_akpar_gk = convol(bracket_akpar_gk)
    total_bracket = Bracket(Dxapar,Dyapar,dx_bracket_akpar_gk + rhos_de*sqrt(ngtot*1.0)*Dxgm,
        dy_bracket_akpar_gk + rhos_de*sqrt(ngtot*1.0)*Dygm)
    
    fglast = -bracket_phik_gk + total_bracket

    return fglast

end

