include("constants.jl")
include("grids.jl") #assuming this will exist

# Need to resolve kperp array from grids

function exp_nu(i, j, niu2, dti)
    exp(-( niu*kperp[i, j]^2 + niu2*kperp[i, j]^(2*hyper_order))*dti)
end

function exp_ng(ng,hyper_nuei,dti)
    exp(-(ng*nu_ei+ng^(2*hyper_morder)*hyper_nuei)*dti)
end

function exp_eta(i, j, res2,dti)
    exp(-(res*kperp[i, j]^2+res2*kperp[i, j]^(2*hyper_order))*dti/(1.0+kperp[i, j]^2*de^2))