#  use Brackets,   only: funcne_i, funcAkpar_i, funcg2, funcgm, func_lastg, bracket_3
include("constants.jl")
print(lambda)
function Funcne_i(dxphi, dyphi, dxne, dyne, dxapar, dyapar, dxuepar, dyuepar, fne, braakparuekpar)
    # assuming the necessary constants are defined in a module called "constants"
    nkx_par = constants.nkx_par
    nly_par = constants.nly_par
    nlz_par = constants.nlz_par
    nlx = constants.nlx
    nly_par = constants.nly_par
    nky = constants.nky

    # I don't think you need to declare them like that in julia ^, once you do include(constants.jl) you can just use the name directly

    braphiknek = Array{ComplexF64}(undef, nky, nkx_par, nlz_par)
    for k in 1:nlz_par, i in 1:nkx_par, j in 1:nky
        braphiknek[j, i, k] = dxphi[i, j, k] + im * dyphi[i, j, k] + dxne[i, j, k] + im * dyne[i, j, k]
    end

    braakparuekpar = Array{ComplexF64}(undef, nky, nkx_par, nlz_par)
    for k in 1:nlz_par, i in 1:nkx_par, j in 1:nky
        braakparuekpar[j, i, k] = dxapar[i, j, k] + im * dyapar[i, j, k] + dxuepar[i, j, k] + im * dyuepar[i, j, k]
    end

    for k in 1:nlz_par, i in 1:nkx_par, j in 1:nky
        fne[j, i, k] = -braphiknek[j, i, k] + braakparuekpar[j, i, k] # assuming the commented line is not required
    end
    return nothing  # assuming the "braakparuekpar" output argument was modified in-place
end

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

function Bracket_3(dxF, dyF, dxG, dyG)
    # this function calculates the "bracket" of two quantities, i.e., 
    # dxk/dx*dyk/dy-dxk/dy*dyk/dx. inputs are xk and yk, the quantities in
    # k space. Returns braxyk, the value of the bracket also in k-space.

    braxyk = zeros(nky, nkx_par, nlz_par)

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


function Funcne_i(Dxphi, Dyphi, Dxne, Dyne, DxApar, DyApar, Dxuepar, Dyuepar)
    fne = zeros(ComplexF64, nky, nkx_par, nlz_par)

    #ya need to figure out bracket!
    braphiknek = bracket(Dxphi, Dyphi, Dxne, Dyne)
    braakparuekpar = bracket(DxApar, DyApar, Dxuepar, Dyuepar)

    @. fne = -braphiknek + braakparuekpar

    return fne, braakparuekpar
end

function FuncAkpar_i(DxApar,DyApar,Dxphi,Dyphi,Dxne,Dyne,Dxuepar,Dyuepar,Dxg2,Dyg2,akpar,t)
    FApar = zeros(ComplexF64, nky, nkx_par, nlz_par)

    tempx = Dxphi - rhos^2*(Dxne +sqrt(2.0)*Dxg2 )
    tempy = Dyphi - rhos^2*(Dyne +sqrt(2.0)*Dyg2 )

    braakparphik = bracket(DxApar, DyApar, tempx, tempy)
    braakparuekpar = bracket(Dxuepar, Dyuepar, Dxphi, Dyphi)

    # these lines have like some random constants that are never defined? or commente dout? idk whats going on there 
    # also braakparuekpar never used in this subroutine. 
    for k = 1:nlz_par
        for i = 1:nkx_par
            for j = 1:nky
                fapar[j, i, k] = 1.0/(1.0+kperp[j, i]^2*de^2)*
                    (braakparphik[j, i, k] - de^2*brauekparphik[j, i, k] -
                     notanj*1.0/sqrt(2.0)*rhos_de*rhoe_LTe*(0.0,1.0)*ky[j]*akpar[j, i, k])
            end
        end
    end


    return FApar

end 

function Funcg2(Dxg2, Dyg2, Dxphi, Dyphi, Dxapar, Dyapar, Dxg3, Dyg3, braakparuekpar, phik)
    nlx, nly_par, nlz_par = size(Dxapar)
    nky, nkx_par, nlz_par = size(braakparuekpar)

    fg2 = zeros(ComplexF64, nky, nkx_par, nlz_par)

    brag2phik = bracket(Dxg2, Dyg2, Dxphi, Dyphi)
    braakparg3 = bracket(Dxapar, Dyapar, Dxg3, Dyg3)

    for k in 1:nlz_par, i in 1:nkx_par, j in 1:nky
        fg2[j,i,k] = brag2phik[j,i,k] + sqrt(gmin+1.0) * rhos_de * braakparg3[j,i,k] +
                     notanj * sqrt(2.0) * braakparuekpar[j,i,k] -
                     notanj * 1.0/(2.0*rhos_de) * rhoe_LTe*(0.0,1.0) * ky[j] * phik[j,i,k]
    end

    return fg2
end

# ok now funcgm is used by two interfaces. not gonna do it now since I think they r related to the dimensionality and maybe we only need to do one. 

function Func_lastg(hyper_nuei, niu2, Dxgm, Dygm, Dxg, Dyg, Dxphi, Dyphi, Dxapar, Dyapar)
    # Calculates g_M using the analytical nonlinear closure
    # As is, only does the perp part
    
    # Use statements
    # Not sure about nlx, nly_par, nkx_par, nlz_par, nky, ngot, nu_ei, kperp, niu, and hyper_morder// should be from constants
    # Define constants and modules if necessary (the grid stuff idk)
    # Omitted for brevity
    
    # Input:
    # hyper_nuei: real
    # niu2: real
    # Dxgm: real[nlx, nly_par, nlz_par]
    # Dygm: real[nlx, nly_par, nlz_par]
    # Dxg: real[nlx, nly_par, nlz_par]
    # Dyg: real[nlx, nly_par, nlz_par]
    # Dxphi: real[nlx, nly_par, nlz_par]
    # Dyphi: real[nlx, nly_par, nlz_par]
    # Dxapar: real[nlx, nly_par, nlz_par]
    # Dyapar: real[nlx, nly_par, nlz_par]
    
    # Output:
    # fglast: complex[nky, nkx_par, nlz_par]

    brafikg = Bracket(Dxphi, Dyphi, Dxg, Dyg)
    braakparg = Bracket(Dxapar, Dyapar, Dxg, Dyg)

    for k = 1:nlz_par
        for i = 1:nkx_par
            for j = 1:nky
                braakparg[j, i, k] = (rhos_de^2 * (ngtot+1)) / ((ngtot+1) * nu_ei +
                    (ngtot+1)^(2*hyper_morder) * hyper_nuei +
                    niu * kperp[j, i]^2 + niu2*kperp[j, i]^(2*hyper_order)) * braakparg[j, i, k]
            end
        end
    end

    #some convol stuff. seems to be an interface in diag.f90
    #call another bracket
end

