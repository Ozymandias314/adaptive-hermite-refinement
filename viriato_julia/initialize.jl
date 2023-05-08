include("constants.jl")
include("grid.jl")
include("transforms.jl")

function init_perturb()
    field = Array{Real}(undef, nlx, nly_par)
    field .= 0
    
    # fieldk = Array{Complex}(undef, nlx, nly_par, nlz_par)
    # fieldk .= 0
    # nmax = nlz/4
    
    if perturb_type=="none"
        # emptiness
    end 

    # not really sure how this works
    # if perturb_type == "allk"
    #     for k in 1:nlz_par
    #         for i in 1:nkx_par
    #             for j in 1:nky
    #                 fieldk[i,j,k] = (1.0,1.0) * cos(2*pi*zz(k)/Lz)
    #             end
    #         end
    #         FFT2d_inv(fieldk(:, :, k), field(:, :, k))
    #     end
    # end
    field = perturb_amp*field
    return field
end

function equilibrium()
    Apar_eq = Array{Real}(undef, nlx, nly)
    phi_eq = similar(Apar_eq)
    AKpar_eq = Array{Complex}(undef, nkx, nky)
    
    phi_eq .= 0.0
    Apar_eq .= 0.0
    AKpar_eq .= 0.0

    if equilib_type=="gaus"
        for i in 1:nlx     
            for j in 1:nly
                Apar_eq[i,j] = a0*exp(-(yy(j)*2*pi*2/ly)^2)*
                                exp(-(xx(i)*2*pi*2/lx)^2)
                phi_eq[i,j]=0
            end do
        end do
        
        #not sure what happens in the FFT
        FFT2d_direct(Apar_eq(:, :), AKpar_eq(:, :)) # This is needed to calc the perturbed spectrum. Not sure this will work if not done in 3D, so be mindful of this. A. Velb 5/26/22
     end

    return Apar_eq, phi_eq