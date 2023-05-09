# ********** CONSTANT DECLARATION *************
# BOX PARAMETERS:
#const three_d = false
#const PI = 3.1415926535 # TODO: Just use Julia Pi?
#commenting all of these out beacuse you read them from the json
#const lx = 1.0  # X box size in units of 2*pi
#const ly = 1.0  # Y box size in units of 2*pi
#const lz = 1.0  # Z box size in units of 2*pi
#const nlx = 16  # number of grid points in X
#const nly = 16  # number of grid points in Y
#const nlz = 16  # number of grid points in Z
#const mm = 2
#const nn = 2
#const g_inc = false  # true included g's in the calculation
const gtot = 0  # number of Hermite moments; set to zero to run without g's
const gmin = 0
#const dim_vec = 2  # dimension of the solution vector. 
             # E.g., if g_inc=false dim_vec=2 (ne and Apar)
#const NPE = 1  # number of processors
#const NPEz = 1
# time parameters:
#const tmax = 1000  # maximum number of iterations	
#const init_aa0_fac = 0.1  # multiplier on the SI operator   
  
#const CFL_frac = 0.25
#const epsilon = 1e-10     # sets the tolerance for the pth iteration at each time step
# const g_epsilon=1.e-3
#const pmax = 1  # maximum number of iterations per time step 
# FLR:      
#const rhoi = 1e-8
#const rhos = 1e-8
#const de = 1e-8
const rhos_de = 0.0
# MHD:
#const small_rhoi = 1e-6  # if rhoi<small_rhoi, code evolves the RMHD eqs.
# ANJOR:
#const anjor = true
#const notanj = 1.0 # TODO: should just set all instances of notanj to 1.0
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
#const rhoe_LTe = 0.0
const kc0 = 0.0
# DIFFUSION:
#const nu_ei = 0.00	  # collisions
#const res = 0.0  # nu_ei*de**2
#const niu = 0.0  # viscosity
#const hyper_coef = 1.0  # 0.25   # set to zero for no hyper-diffusion
#const hyper_coef_g = 0.0  # 0.25   # set to zero for no hyper-diffusion on the g's
#const hyperm_coef = 2.0  # 0.25
#const hyper_fixed = false  # Fixes the hyper coeffs, if false, scales with d
#const hyper_order = 3  # this means it's k^(2*hyper_order)
#const hyper_order_g = 3  # this means it's k^(2*hyper_order_g) in the g eqs
#const hyper_morder = 3  # this means it's k^(2*hyper_order)
#const hyper_eta = 1.0  # 3.0e-8  
#const hyper_nu::Float64 = 1.0 # 3.0e-8
#const hyper_nu_g::Float64 = 0.0 # 0.25
#const hyper_colls_fixed::Bool = false
#const hyper_colls::Float64 = 0.0

#********** EQUILIBRIUM VALUES **************
#const A0::Float64 = 0.0
#const PHI0::Float64 = 0.0 # not used
#const Leq::Float64 = 0.3
# const perturb_A::Float64 = 0.0 # 1.e-5*A0
# const perturb_PHI::Float64 = 0.0 # 1.0e-8 or 1.e-5*PHI0
# const g_perturb::Float64 = 0.0
#const perturb_amp::Float64 = 0.0
#const equilib_type::String = "none" # AVK: Changed default equil to none
#const perturb_type::String = "none"
#*********************************
#**** TURBULENT SOURCE PARAMETERS ****** TODO: Not needed!
# const turb::Bool = true
# # const gturb::Bool = false
# const kfp1::Float64 = 1.0
# const kfp2::Float64 = 2.0
# const kfz1::Float64 = 1.0
# const kfz2::Float64 = 1.0
# const feps::Float64 = 0.0


#Skipped a bunch of stuff related to saving data

#definitions I guess

# TODO:give real definitions! Some things here not needed, others are depended on constants defined above, like nkx,nky, dx,dy.
# vincent: nlx and stuff are above, they used to be upper case whoops

# See Viriato code for clarification? Idk
#nkx::Int64 = nlx/2+1
#nky::Int64 = nly
#nkz::Int64
#x_loc::Int64
#y_loc::Int64
#z_loc::Int64
#k_max::Int64 = nkx - div(nkx, 2)
#kperpmax::Float64 = sqrt(nkx^2 + (nky/2)^2) # confused! changed to float64
#nkx_par::Int64 # TODO: make sure to change all uses of "_par" to just the nkx as this has to do with the parallelization!
#nly_par::Int64
#nlypar_old::Int64
#nlz_par::Int64
#dx::Float64
#dy::Float64
#dz::Float64
#etaz::Float64 = 0.0
#etaz_g::Float64 = 0.0
#scale::Float64

#read data from input file

using JSON

# Load the JSON file
const_data = JSON.parsefile("viriato_julia/example_inputs.json")

# Extract values from the "box_parameters" object
lx = const_data["box_parameters"]["lx"]
ly = const_data["box_parameters"]["ly"]
lz = const_data["box_parameters"]["lz"]
nlx = const_data["box_parameters"]["nlx"]
nly = const_data["box_parameters"]["nly"]
nlz = const_data["box_parameters"]["nlz"]
g_inc = const_data["box_parameters"]["g_inc"]
ngtot = const_data["box_parameters"]["ngtot"]
npe = const_data["box_parameters"]["npe"]
npez = const_data["box_parameters"]["npez"]
three_D = const_data["box_parameters"]["three_D"]

# Extract values from the "time_parameters" object
tmax = const_data["time_parameters"]["tmax"]
init_aa0_fac = const_data["time_parameters"]["init_aa0_fac"]
cfl_frac = const_data["time_parameters"]["cfl_frac"]
epsilon = const_data["time_parameters"]["epsilon"]
pmax = const_data["time_parameters"]["pmax"]

# Extract values from the "flr" object
rhoi = const_data["flr"]["rhoi"]
rhos = const_data["flr"]["rhos"]
de = const_data["flr"]["de"]

# Extract values from the "mhd" object
small_rhoi = const_data["mhd"]["small_rhoi"]

# Extract values from the "diffusion" object
nu_ei = const_data["diffusion"]["nu_ei"]
res = const_data["diffusion"]["res"]
niu = const_data["diffusion"]["niu"]
hyper_coef = const_data["diffusion"]["hyper_coef"]
hyper_coef_g = const_data["diffusion"]["hyper_coef_g"]
hyperm_coef = const_data["diffusion"]["hyperm_coef"]
hyper_fixed = const_data["diffusion"]["hyper_fixed"]
hyper_order = const_data["diffusion"]["hyper_order"]
hyper_order_g = const_data["diffusion"]["hyper_order_g"]
hyper_morder = const_data["diffusion"]["hyper_morder"]
hyper_eta = const_data["diffusion"]["hyper_eta"]
hyper_nu = const_data["diffusion"]["hyper_nu"]
hyper_nu_g = const_data["diffusion"]["hyper_nu_g"]
hyper_colls_fixed = const_data["diffusion"]["hyper_colls_fixed"]
hyper_colls = const_data["diffusion"]["hyper_colls"]

# Extract values from the "equil" object
a0 = const_data["equil"]["a0"]
phi0 = const_data["equil"]["phi0"]
leq = const_data["equil"]["leq"]
perturb_amp = const_data["equil"]["perturb_amp"]
equilib_type = const_data["equil"]["equilib_type"]
perturb_type = const_data["equil"]["perturb_type"]

# Extract values from the "data_sav" object
save_energyfiles = const_data["data_sav"]["save_energyfiles"]
save_checkpoints = const_data["data_sav"]["save_checkpoints"]
save_datafiles = const_data["data_sav"]["save_datafiles"]
restart = const_data["data_sav"]["restart"]
rs_time = const_data["data_sav"]["rs_time"]



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