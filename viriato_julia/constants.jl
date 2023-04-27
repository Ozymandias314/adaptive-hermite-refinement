# ********** CONSTANT DECLARATION *************
# BOX PARAMETERS:
const three_d = true
const PI = 3.1415926535
const Lx = 1.0  # X box size in units of 2*pi
const Ly = 1.0  # Y box size in units of 2*pi
const Lz = 1.0  # Z box size in units of 2*pi
const NLx = 16  # number of grid points in X
const NLy = 16  # number of grid points in Y
const NLz = 16  # number of grid points in Z
const mm = 2
const nn = 2
const g_inc = false  # true included g's in the calculation
const gtot = 0  # number of Hermite moments; set to zero to run without g's
const gmin = 0
const dim_vec = 2  # dimension of the solution vector. 
             # E.g., if g_inc=false dim_vec=2 (ne and Apar)
const NPE = 1  # number of processors
const NPEz = 1
# time parameters:
const tmax = 1000  # maximum number of iterations	
const init_aa0_fac = 0.1  # multiplier on the SI operator   
  
const CFL_frac = 0.25
const epsilon = 1e-10     # sets the tolerance for the pth iteration at each time step
# const g_epsilon=1.e-3
const pmax = 1  # maximum number of iterations per time step 
# FLR:      
const rhoi = 1e-8
const rhos = 1e-8
const de = 1e-8
const rhos_de = 0.0
# MHD:
const small_rhoi = 1e-6  # if rhoi<small_rhoi, code evolves the RMHD eqs.
# ANJOR:
const anjor = true
const notanj = 1.0
const zcharge = 1.0
const tite = 1.0
const betai = 1.0
const plusorminus = 1
const inkx = 1
const inky = 1
const inkz = 1

const lambda = 1.0
const sigma = 1.0
const rhos_diag = 0.0
# Background electron temp gradient:
const rhoe_LTe = 0.0
const kc0 = 0.0
# DIFFUSION:
const nu_ei = 0.00	  # collisions
const res = 0.0  # nu_ei*de**2
const niu = 0.0  # viscosity
const hyper_coef = 1.0  # 0.25   # set to zero for no hyper-diffusion
const hyper_coef_g = 0.0  # 0.25   # set to zero for no hyper-diffusion on the g's
const hyperm_coef = 2.0  # 0.25
const hyper_fixed = false  # Fixes the hyper coeffs, if false, scales with d
const hyper_order = 3  # this means it's k^(2*hyper_order)
const hyper_order_g = 3  # this means it's k^(2*hyper_order_g) in the g eqs
const hyper_morder = 3  # this means it's k^(2*hyper_order)
const hyper_eta = 1.0  # 3.0e-8  
const hyper_nu::Float64 = 1.0 # 3.0e-8
const hyper_nu_g::Float64 = 0.0 # 0.25
const hyper_colls_fixed::Bool = false
const hyper_colls::Float64 = 0.0

#********** EQUILIBRIUM VALUES **************
const A0::Float64 = 0.0
const PHI0::Float64 = 0.0 # not used
const Leq::Float64 = 0.3
# const perturb_A::Float64 = 0.0 # 1.e-5*A0
# const perturb_PHI::Float64 = 0.0 # 1.0e-8 or 1.e-5*PHI0
# const g_perturb::Float64 = 0.0
const perturb_amp::Float64 = 0.0
const equilib_type::String = "none" # AVK: Changed default equil to none
const perturb_type::String = "none"
#*********************************
#**** TURBULENT SOURCE PARAMETERS ******
const turb::Bool = true
# const gturb::Bool = false
const kfp1::Float64 = 1.0
const kfp2::Float64 = 2.0
const kfz1::Float64 = 1.0
const kfz2::Float64 = 1.0
const feps::Float64 = 0.0


#Skipped a bunch of stuff related to saving data

#definitions I guess

nkx::Int64
nky::Int64
nkz::Int64
x_loc::Int64
y_loc::Int64
z_loc::Int64
k_max::Int64
kperpmax::Int64
nkx_par::Int64
nly_par::Int64
nlypar_old::Int64
nlz_par::Int64
dx::Float64
dy::Float64
dz::Float64
etaz::Float64 = 0.0
etaz_g::Float64 = 0.0
scale::Float64

# Definitions
nkx = div(nlx, 2) + 1
nky = nly
nkz = 2 * div(nlz - 1, 3) + 1
x_loc = div(nlx, 2) + 1
y_loc = div(nly, 2) + 1
z_loc = div(nlz, 2) + 1
k_max = nkx - div(nkx, 2)  # maximum k after filtering
kperpmax = sqrt(nkx^2 + div(nky, 2)^2)
nkx_par = div(nkx - 1, npe) + 1
nly_par = div(nly - 1, npe) + 1
nlz_par = div(nlz - 1, npez) + 1
dx = lx / nkx
dy = ly / nky
dz = lz / nlz
etaz = etaz * dz^2  # AVK z viscosity will act on kz .ge. nlz/3 == dealiasing?
etaz_g = etaz_g * dz^2  # AVK z viscosity will act on kz .ge. nlz/3 == dealiasing?
scale = 1.0 / (nlx * nly)  # scale factor for FFTs

j1 = div(kperp0 * ly, lx) + 1
j2 = nky - j1 + 2