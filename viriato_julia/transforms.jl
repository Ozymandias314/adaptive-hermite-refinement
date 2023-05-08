using LinearAlgebra,FFTW

#need to make plans global?

# 2D FFT
function FFT2d_direct(array::Array{Float64},first_fft::Bool)
    
    if first_fft
        rcfft_plan = plan_rfft(array)
        first_fft = false
    end
    
    array_k = rcfft_plan*array
end

# 2D IFFT
function FFT2d_inv(array_k::Array{ComplexF64},first_ifft::Bool)
    

    if first_ifft
        irfft_plan = plan_irfft(array_k,nlx) 
        first_ifft = false
    end

    array = irfft_plan*array_k
end