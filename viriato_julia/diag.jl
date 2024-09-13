include("constants.jl")
include("grid.jl")
using FFTW


function print_cpp(arr::Matrix{ComplexF64})
    for i in 1:size(arr, 1)
        row_str = ""
        for j in 1:size(arr, 2)
            real_part = round(real(arr[i, j]), digits=6)
            imag_part = round(imag(arr[i, j]), digits=6)
            row_str *= "($(real_part), $(imag_part)) "
        end
        println(row_str)
    end
end

function print_cpp(arr::Matrix{Float64})
    for i in 1:size(arr, 1)
        row_str = ""
        for j in 1:size(arr, 2)
            real_part = sprintf1("%.5e",arr[i, j])
            row_str *= "$(real_part) "
        end
        println(row_str)
    end
end


function convol(Fk)
    #print(Fk[32,32],"\n")
    DxF = Array{Float64}(undef, nlx, nly)
    DyF = similar(DxF)

    Fk_ikx = similar(Fk)
    Fk_iky = similar(Fk)

    for i in 1:nkx
        for j in 1:nky
            Fk_ikx[i,j] = im*kx(i)*Fk[i,j] / nlx / nly
            Fk_iky[i,j] = im*ky(j)*Fk[i,j] / nlx / nly
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