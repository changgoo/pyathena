<comment>
problem         = tigress_default # 
author          = C.-G. Kim       # 
journal         = TIGRESS (Kim & Ostriker 2017) # 
config          =                 # 

<job>
problem_id      = R2              # problem ID: basename of output filenames
maxout          = 5               # Output blocks number from 1 -> maxout
num_domains     = 1               # number of Domains in Mesh

<log>
file_open       = 0               # 
iflush          = 1               # Default Value
out_level       = 0               # Default Value
err_level       = 0               # Default Value

<output1>
out_fmt         = hst             # History data dump
dt              = 0.01            # time increment between outputs
time            = 0.0             # Next Output Time
num             = 0               # Next Output Number
level           = -1              # Default Value
domain          = -1              # Default Value
id              = out1            # Default Value
out             = cons            # Default Value

<output2>
out_fmt         = vtk             # Binary data dump
out             = prim            # 
dt              = 0.2             # time increment between outputs
time            = 0.0             # Next Output Time
num             = 0               # Next Output Number
level           = -1              # Default Value
domain          = -1              # Default Value
id              = out2            # Default Value

<output3>
out_fmt         = starpar_vtk     # 
dt              = 0.2             # 
time            = 0.0             # Next Output Time
num             = 0               # Next Output Number
level           = -1              # Default Value
domain          = -1              # Default Value
id              = out3            # Default Value
out             = cons            # Default Value

<output4>
out_fmt         = rst             # Restart data dump
dt              = 11.904762       # time increment between outputs
time            = 0.0             # Next Output Time
num             = 0               # Next Output Number
level           = -1              # Default Value
domain          = -1              # Default Value
id              = out4            # Default Value
out             = cons            # Default Value

<output5>
out_fmt         = zprof           # 
dt              = 0.2             # time increment between outputs

<time>
grav_no         = 1               # 
cour_no         = 0.3             # The Courant, Friedrichs, & Lewy (CFL) Number
nlim            = -1              # cycle limit
tlim            = 100.            # time limit
time            = 0.0             # Current Simulation Time
nstep           = 0               # Current Simulation Time Step

<domain1>
level           = 0               # refinement level this Domain (root=0)
Nx1             = 128             # Number of zones in X-direction
x1min           = -128            # minimum value of X
x1max           = 128             # maximum value of X
bc_ix1          = 4               # inner-I (X) boundary flag
bc_ox1          = 4               # outer-I (X) boundary flag
NGrid_x1        = 8               # 
Nx2             = 128             # Number of zones in Y-direction
x2min           = -128            # minimum value of Y
x2max           = 128             # maximum value of Y
bc_ix2          = 4               # inner-J (Y) boundary flag
bc_ox2          = 4               # outer-J (Y) boundary flag
NGrid_x2        = 8               # 
Nx3             = 1792            # Number of zones in X3-direction
x3min           = -1792           # minimum value of X3
x3max           = 1792            # maximum value of X3
bc_ix3          = 2               # inner-K (Z) boundary flag
bc_ox3          = 2               # outer-K (Z) boundary flag
NGrid_x3        = 28              # 
AutoWithNProc   = 0               # 

<problem>
gamma           = 1.66666667      # gamma = C_p/C_v
amp             = 0.0             # drho/rho <= amp
kappa           = 1.e7            # conductivity in unit of c.g.s.
iprob           = 5               # 1: eigen mode/2: random perturbation/3: power spectrum
surf            = 144.709         # gas surface density
sz0             = 30              # initial velocity dispersion
n0              = 14.44           # initial midplane density (sets real gas surface density)
beta            = 10              # P_th/P_mag
ispect          = 1               # 1: power law/2: G&O form
expo            = 2.000000        # exponent of power law. for ispect=1
kpeak           = 4               # peak wave number. for ispect=2
klow            = 1               # cutoff wavenumber
khigh           = 128             # cutoff wavenumber
Omega           = 0.0908776       # angular speed of Galactic rotation (km/s/kpc)
qshear          = 1.0             # shear parameter -dlnOmega/dlnR
SurfS           = 462.977         # stellar surface density for KG89 model (Msun/pc^2)
zstar           = 245             # stellar disk scale height for KG89 model (pc)
rhodm           = 0.0791234       # dark matter density at the midplane for KG89 model (Msun/pc^3)
R0              = 2000            # Galactocentric radius (pc)
dedt            = 1.00            # Amplitude of driven turbulence: Edot per area (Lsun/pc^2)
dtdrive         = 0.1             # Time interval for new perturbation (Myr)
driving_time    = 20.0            # Shut off time for driven turbulence (Myr)
heat_ratio0     = 50.0            # floor for the heating up to driving_time
dt_init_cluster = 10.0            # age range of initial clusters
H_init_cluster  = 10.0            # scale height of initial clusters
randomSN_rate   = 0.0             # rate of the random SN component
randomSN_height = 100.0             # scale height of the random SN component

<feedback>
SFUV0           = 42.6351         # G/G0=SFUV/SFUV0
fbin            = 0.6666667       # fraction of binary runaways
fdyn            = 0.0             # fraction of dynamical runaways
Mej             = 10.0            # mass of ejecta (+ pre winds)
eps             = 1.0             # m_star/m_sp
iSN             = 5               # mode of SN: 1-th+mom; 2-th+ST+mom; 3-th; 4-mom
Ia_amp		= 0.0		  # Amplitude of type Ia SNe
