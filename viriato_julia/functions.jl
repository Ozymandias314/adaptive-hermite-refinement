include("constants.jl")
include("grids.jl") #assuming this will exist

# Need to resolve kperp array from grids

function exp_nu(j, i, niu2, dti)
    exp(-( niu*kperp[j,i]^2 + niu2*kperp[j,i]^(2*hyper_order))*dti)
end

function exp_ng(dti)
    exp(-(res*kperp[j,i]^2+res2*kperp[j,i]^(2*hyper_order))*dti/(1.0+kperp[j,i]^2*de^2))
end

function exp_eta(res2,dti)
    exp(-(res*kperp[j,i]^2+res2*kperp[j,i]^(2*hyper_order))*dti/(1.0+kperp[j,i]^2*de^2))