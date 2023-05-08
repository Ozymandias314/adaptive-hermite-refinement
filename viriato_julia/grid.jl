include("constants.jl") # Need to include constants everywhere?


function kperp(i::Int,j::Int)
    kperp = sqrt(ky(j)^2+kx(i)^2)
end

#not strictly needed, but maybe later for parallelization?
function kx(i::Int)
    kx = i
end

# Get ky to respect the order of frequencies returned by fft WITHOUT any FFTSHIFT calls
function ky(j::Int)
    if i <= nky/2+1
        ky = (j-1)*lx/ly
    else
        ky = (j-nky-1)*lx/ly
    end
end

function Γ₀(x::Float64)
    ax = abs(x)
    ans = zero(x)
    y = zero(x)
    
    if ax < 3.75
        y = (x / 3.75)^2
        ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 +
            y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))))
    else
        y = 3.75 / ax
        ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y * (0.1328592e-1 +
            y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2 +
            y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 +
            y * 0.392377e-2))))))))
    end
    
end