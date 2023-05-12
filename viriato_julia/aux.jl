include("constants.jl")
include("grid.jl")

function omegakaw(bperp_max::Real)
    kperp_dum::Real = 0.0
    first::Bool = true # TODO: This isnt quite right. In Fortran, have the save statement, which saves this bool for future calls of this function. Right now there is unecessary computation

    if first
        kperp_dum = sqrt(ky(nky÷2+1)^2 + (nkx*1.0)^2)
        first = false
    end

    if rhoi < small_rhoi
        omega_kaw = sqrt(1.0 + kperp_dum^2*(3.0/4.0*rhoi^2 + rhos^2)) * ky(nky÷2+1) * bperp_max / (1.0 + kperp_dum^2 * de^2)
        omega_kaw_nl = sqrt(1.0 + kperp_dum^2*(3.0/4.0*rhoi^2 + rhos^2)) * kperp_dum * bperp_max / (1.0 + kperp_dum^2 * de^2)
    else
        omega_kaw = kperp_dum * sqrt(rhos^2 - rhoi^2 / (Γ₀(0.5*kperp_dum^2*rhoi^2)-1.)) *
            ky(nky÷2+1) * bperp_max / sqrt(1.0 + kperp_dum^2 * de^2)
        omega_kaw_nl = kperp_dum^2 * sqrt(rhos^2 - rhoi^2 / (Γ₀(0.5*kperp_dum^2*rhoi^2)-1.)) *
            bperp_max / (1.0 + kperp_dum^2 * de^2)
    end

    return omega_kaw
end

#Note: @inbounds macro is used to eliminate bounds checking inside the loops for performance optimization.

function phi_pot(nek::Array{ComplexF64,2})
    @inbounds begin
        phiK = zeros(ComplexF64,nkx,nky)

        if (rhoi < small_rhoi)
            for i in 1:nkx
                for j in 1:nky
                    phiK[i,j] = -nek[i,j]/(kperp(i,j)^2)
                end
            end
        else
            for i in 1:nkx
                for j in 1:nky
                    phiK[i,j] = rhoi^2*0.5/(Γ₀(kperp(i,j)^2*rhoi^2*0.5)-1.0)*nek[i,j]
                end
            end
        end
    end
    return phiK
end

function func_semi_implicit_operator(dti::Real, bperp_max::Real, aa0::Real)
    SI_oper = Array{ComplexF64}(undef, nkx, nky)
    if rhoi <= small_rhoi
        for i in 1:nkx
            for j in 1:nky
                SI_oper[i,j] = aa0^2*(1+kperp(i,j)^2*(3/4*rhoi^2+rhos^2))*kperp(i,j)^2*bperp_max^2*
                    dti^2/(1.0+kperp(i,j)^2*de^2) 
            end
        end
    else
        for i in 1:nkx
            for j in 1:nky
                SI_oper[i,j] = aa0^2*(3*rhos^2-rhoi^2/(Γ₀(0.5*kperp(i,j)^2*rhoi^2)-1))*
                                kperp(i,j)^4*bperp_max^2*
                                dti^2/(1.0+kperp(i,j)^2*de^2) 
            end
        end
    end
    return SI_oper
end

function dtnext(relative_error::Float64,dti_temp::Float64,noinc::Bool,dti::Float64)
    if noinc
        inc_fac = 1.0
        noinc = false # need to make noinc global
    else
        inc_fac = 1.08
    end

    if relative_error < 0.8*epsilon
        if dti_temp < inc_fac*dti
            dti = dti_temp
        else
            dti = inc_fac*dti
        end
    else
        dti = min(dti_temp,dti)
    end

    return dti,noinc
end