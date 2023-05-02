include("constants.jl")
using FFTW

function Convol(Fk)
    nkx_par = constants.nkx_par
    nky = constants.nky
    nlx = constants.nlx
    nly_par = constants.nly_par
    nlz_par = constants.nlz_par
    
    # Does Fk need to be defined again?
    # Fk = Array{ComplexF64}(undef, nky, nkx_par, nlz_par)
    DxF = Array{Float64}(undef, nlx, nly_par, nlz_par)
    DyF = similar(DxF)

    Fk_ikx = similar(Fk)
    Fk_iky = similar(Fk)

    for k in 1:nlz_par
        for i in 1:nkx_par
            for j in 1:nky
                Fk_ikx[j,i,k] = im*kx[i]*Fk[j,i,k]
                Fk_iky[j,i,k] = im*ky[i]*Fk[j,i,k]
            end
        end
    end

    # Assuming DxF, DyF are iFFT of Fk_ikx
    DxF = ifft(Fk_ikx)
    DyF = ifft(Fk_iky)

    return DxF, DyF
end

function flows(dxfi,dyfi,dxne,dyne)
    # Assuming you don't have to define inputs again
    vex = Array{Float64}(undef, nlx, nly_par, nlz_par)
    vey = similar(vex)
    vrhosx = similar(vex)
    vrhosy = similar(vex)

    if rhoi < small_rhoi
        vex = -dyfi
        vey = dxfi
        vrhosx = 0.0
        vrhosy = 0.0
    else
        vex = -dyfi
        vey = dxfi
        vrhosx = -rhos^2*dxne
        vrhosy = rhos^2*dyne
    end
    return vex, vey, vrhosx, vrhosy
end

function b_field(dxapar,dyapar)
    bx = Array{Float64}(undef, nlx, nly_par, nlz_par)
    by = similar(bx)
    bx = dyapar
    by = -dxapar
    return bx, by