using LinearAlgebra,FFTW

#need to make plans global?
# make plans here. They are resused whenever we call FFT, ie the FFTs are always the same size. These will be done the first time the file is included, at the beginning of REGK. 
# Note that nkx = nlx/2 +1, nky = nly, since we do the real to complex fft and therefore only need one half of the FFT in the x direction.
array = Array{Float64}(undef,nlx,nly)
array_k = Array{ComplexF64}(undef,nkx,nky) # 
rcfft_plan = plan_rfft(array)
brfft_plan = plan_brfft(array_k,nlx) # Supply the first dimension of the real space array that was shrunk

# 2D FFT
function FFT2d_direct(array::Array{Float64})
    array_k = rcfft_plan*array 
    # Do Fourier filtering
    if hou_li_filter
        for i in 1:nkx
            array_k[i,:] *= exp(-36*(kx(i)/nkx)^36)
        end
        for j in 1:nky
            array_k[:,j] *= exp(-36*(ky(j)/nky)^36)
        end
    else 
        for i in 1:nkx
            if kx(i)/nkx > 2.0/3.0
                array_k[i,:] .= 0.0
            end
            
        end
            
        for j in 1:nky
            if ky(j)/nky > 2.0/3.0
                array_k[:,j] .= 0.0
            end
        end
    end
    return array_k
end

# 2D IFFT
function FFT2d_inv(array_k::Array{ComplexF64})
    array = brfft_plan*array_k
end