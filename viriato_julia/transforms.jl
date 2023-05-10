using LinearAlgebra,FFTW

#need to make plans global?
# make plans here. They are resused whenever we call FFT, ie the FFTs are always the same size. These will be done the first time the file is included, at the beginning of REGK. 
# Note that nkx = nlx/2 +1, nky = nly, since we do the real to complex fft and therefore only need one half of the FFT in the x direction.
array = Array{Float64}(undef,nlx,nly)
array_k = Array{ComplexF64}(undef,nkx,nky) # 
rcfft_plan = plan_rfft(array)
irfft_plan = plan_irfft(array_k,nlx)

# 2D FFT
function FFT2d_direct(array::Array{Float64})
    array_k = rcfft_plan*array/sqrt(nlx*nly) # Must normalize by number of modes! Luka made a comment that we might want to do this in the bracket funciton instead of here, but can change this later. 
    # Do Fourier filtering
    if hou_li_filter
        for i in 1:nkx
            for j in 1:nky
                array_k[i,j] *= exp(-36*(kperp(i,j)/kperpmax)^36)
            end
        end
    else 
        for i in 1:nkx
            for j in 1:nky
                if kperp(i,j)/kperpmax > 2.0/3.0
                    array_k[i,j] *= 0.0
                end 
            end
        end

    end
    return array_k
end

# 2D IFFT
function FFT2d_inv(array_k::Array{ComplexF64})
    array = irfft_plan*array_k/sqrt(nlx*nly)
end