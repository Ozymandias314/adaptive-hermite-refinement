include("constants.jl")
include("grid.jl")
include("transforms.jl")

function init_perturb()
    field = Array{Float64}(undef, nlx, nly)
    field .= 0
    
    fieldk = Array{ComplexF64}(undef, nkx, nky)
    fieldk .= 0
    # nmax = nlz/4
    
    if perturb_type=="none"
        # emptiness
    end 

    # not really sure how this works
    if perturb_type == "allk"
        fieldk .= 1
        field = FFT2d_inv(fieldk)
    end
    field = perturb_amp*field
    return field
end

function equilibrium()
    Apar_eq = Array{Float64}(undef, nlx, nly)
    phi_eq = similar(Apar_eq)
    phi_eq .= 0.0
    Apar_eq .= 0.0
    if equilib_type=="gaus"
        for i in 1:nlx     
            for j in 1:nly
                Apar_eq[i,j] = a0*exp(-(yy(j)*2*pi*2/ly)^2)*
                                exp(-(xx(i)*2*pi*2/lx)^2)
                phi_eq[i,j]=0
            end 
        end
    end
    if equilib_type=="OT01"
        for i in 1:nlx
            for j in 1:nly
                Apar_eq[i,j]=cos(4*pi*xx(i)/lx)+2*cos(2*pi*yy(j)/ly)
                phi_eq[i,j]=-2*(cos(2*pi*xx(i)/lx)+cos(2*pi*yy(j)/ly))
            end
        end
    end
    if equilib_type=="tear"
        for i in 1:nlx
            for j in 1:nly
                Apar_eq[i,j] = -a0/(cosh(xx(i))^2)*1/(2*tanh(pi)^2-tanh(2*pi))*(tanh(xx(i)-pi)^2+tanh(xx(i)+pi)^2-tanh(2*pi)^2)
                phi_eq .= 0
            end
        end
    end
    return Apar_eq, phi_eq
end