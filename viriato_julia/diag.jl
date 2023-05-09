include("constants.jl")
include("grid.jl")
using FFTW

function convol(Fk)
    DxF = Array{Float64}(undef, nlx, nly)
    DyF = similar(DxF)

    Fk_ikx = similar(Fk)
    Fk_iky = similar(Fk)

    for i in 1:nkx
        for j in 1:nky
            Fk_ikx[i,j] = im*kx(i)*Fk[i,j]
            Fk_iky[i,j] = im*ky(j)*Fk[i,j]
        end
    end

    # Assuming DxF, DyF are iFFT of Fk_ikx
    DxF = FFT2d_inv(Fk_ikx)
    DyF = FFT2d_inv(Fk_iky)

    return DxF, DyF
end

function flows(dxfi,dyfi,dxne,dyne)
    vex = Array{Float64}(undef, nlx, nly)
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

function bfield(dxapar,dyapar)
    bx = Array{Float64}(undef, nlx, nly)
    by = similar(bx)
    bx = dyapar
    by = -dxapar
    return bx, by
end