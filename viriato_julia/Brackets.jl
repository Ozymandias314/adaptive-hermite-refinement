#  use Brackets,   only: funcne_i, funcAkpar_i, funcg2, funcgm, func_lastg, bracket_3

function Funcne_i(dxphi, dyphi, dxne, dyne, dxapar, dyapar, dxuepar, dyuepar, fne, braakparuekpar)
    using Main.constants  # assuming the necessary constants are defined in a module called "constants"
    nkx_par = constants.nkx_par
    nly_par = constants.nly_par
    nlz_par = constants.nlz_par
    nlx = constants.nlx
    nly_par = constants.nly_par
    nky = constants.nky

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