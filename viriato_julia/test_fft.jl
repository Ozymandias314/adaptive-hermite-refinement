using LinearAlgebra, FFTW, Plots
include("constants.jl")
include("grid.jl")

test = zeros(Float64,nlx)
d_testk = zeros(ComplexF64,nkx)
for i =1:nlx 
    test[i] = sin(xx(i))
end

testk = rfft(test)
for i = 1:nkx
    d_testk[i] = im*kx(i)*testk[i]
end

d_testk_inv = irfft(d_testk,nlx)
println(xx(1))
#plot(test)
plot(d_testk_inv)
#plot(abs.(testk))
#println("here")
#plot(abs.(d_testk))