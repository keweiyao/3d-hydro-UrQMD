import numpy as np
import sys

def generate_ic_input(	event_index, 
	        	proj,
        		tarj,
        		mult_fluct,
        		nuclen_wdith,
       		 	pp_cross_section,
        		entropy_normalization,
        		rapidity_mean_coeff,
        		rapidity_width,
        		rapidity_skew_coeff,
        		rapidity_kurtosis,
        		Jacobian,
        		ic_x_max,
        		ic_x_step,
        		ix_eta_max,
        		ic_eta_N
		):
	outstring = """
# parameter set number = %d
# specify the projectile option twice
projectile = %s
projectile = %s
number-events = 1

# don't print event properties to stdout, save to HDF5
quiet = true
output = ic.hdf5

reduced-thickness = 0.0

# 1.4 at LHC; > 3.0 at RHIC (The larger this parameter, the smaller the fluctuation)
fluctuation = %f

# 0.43 at LHC; ~0.5 at RHIC
nucleon-width = %f

# 3.98 at 200 GeV; 6.4 at 2.76 TeV; 7.0 at 5.02 TeV
cross-section = %f

normalization = %f

rapidity-mean = %f
rapidity-width = %f
rapidity-skew = %f
rapidity-kurtosis = %f

<pt2/mt2> = %f

grid-max = %f
grid-step = %f
eta-max = %f
eta-step = %i
switch-3d = 1
out-3d = 1
"""%(	event_index,
	proj,
	tarj,
	mult_fluct, 
	nuclen_wdith,
	pp_cross_section,
	entropy_normalization,
	rapidity_mean_coeff,
	rapidity_width,
	rapidity_skew_coeff,
	rapidity_kurtosis,
	Jacobian,
	ic_x_max,
	ic_x_step,
	ix_eta_max,
	ic_eta_N
	)	
	outfile = open("trento-config.txt",'w')
	outfile.write(outstring)
	outfile.close()


def generate_hydro_input(
	out_dir,
	eta_s_min,
	eta_s_QGP_slope,
	eta_s_HRG,
	zeta_s,
	e_critical,
	entropy_scale,
	hydro_nx,
	hydro_ny,
	hydro_nz,
	hydro_xmin,
	hydro_xman,
	hydro_ymin,
	hydro_ymax,
	hydro_zmin,
	hydro_zmax,
	tau_0,
	tau_max,
	dtau
     ):
	outstring  = """
outputDir       %s

etaS_min         %f       ! eta/s value
etaS_slope_QGP   %f       ! QGP slope
etaS_slope_HRG   %f        ! HRG slope

zetaS           %f

e_crit  %f    ! criterion for surface finding
epsilon0       %f      ! Norm factor, 96.6 @ 2.76TeV, 128 @ 5.02 TeV

nx           %i            ! number of cells in X direction
ny           %i            ! number of cells in Y direction
nz           %i             ! number of cells in eta direction
xmin         %f            ! coordinate of the first cell
xmax         %f            ! coordinate of the last cell
ymin         %f
ymax         %f
etamin       %f
etamax       %f

icInputFile     ic.hdf5
event_number    0
tau0    %f     ! starting proper time
tauMax       %f           ! proper time to stop hydro
dtau         %f           ! timestep
"""%(out_dir,
     eta_s_min,
     eta_s_QGP_slope,
     eta_s_HRG,
     zeta_s,
     e_critical,
     entropy_scale,
     hydro_nx,
     hydro_ny,
     hydro_nz,
     hydro_xmin,
     hydro_xman,
     hydro_ymin,
     hydro_ymax,
     hydro_zmin,
     hydro_zmax,
     tau_0,
     tau_max,
     dtau
     )
	outfile = open("vhlle-config.txt",'w')
  	outfile.write(outstring)
     	outfile.close()

def generate_escape_input(output_dir, number_of_samples, e_critical):
	outstring = """
surface freezeout.dat
spectra_dir     %s
number_of_events        %i
rescatter        1
weakContribution 0
shear            1
ecrit   %f
"""%(output_dir, number_of_samples, e_critical)
	outfile = open("escape-config.txt",'w')
	outfile.write(outstring)
	outfile.close()

def generate_stat_input(   m_eta_L,
        m_eta_H, 
        m_pt_L, 
        m_pt_H, 
        dn_eta_L, 
        dn_eta_H, 
        dn_eta_bins, 
        dn_pt_L, 
        dn_pt_H, 
        q_eta_L, 
        q_eta_H,
        q_eta_bins, 
        q_pt_L, 
        q_pt_H):
	outstring = """
inputfile = plist.dat
outputfile = stat.dat

mult = 1
m_eta_L = %f
m_eta_H = %f
m_pt_L = %f
m_pt_H = %f

dndy = 1
dn_eta_L = %f
dn_eta_H = %f
dn_eta_bins = %d
dn_pt_L = %f
dn_pt_H = %f

qns = 1
q_eta_L = %f
q_eta_H = %f
q_eta_bins = %d
q_pt_L = %f
q_pt_H = %f
"""%(	m_eta_L, 
	m_eta_H, 
	m_pt_L, 
	m_pt_H,
	dn_eta_L, 
	dn_eta_H, 
	dn_eta_bins, 
	dn_pt_L, 
	dn_pt_H, 
	q_eta_L, 
	q_eta_H,
	q_eta_bins, 
	q_pt_L, 
	q_pt_H)
	outfile = open("stat-config.txt",'w')
        outfile.write(outstring)
        outfile.close()

if __name__ == "__main__":
	print "generate trentro IC config file"
	generate_ic_input(0, "Pb", "Pb", 2.8, 0.45, 6.4, 1.0, 0.75, 2.8, 0.2, 0.3, 0.65, 10.0, 0.2, 15.0, 64)
#	generate_ic_input(0, "p", "Pb", 3.2, 0.5, 7.0, 1.0, 0.75, 3.0, 0.2, 0.3, 0.65, 5.0, 0.1, 15.0, 64)
	print "generate vhlle config file"
	generate_hydro_input("vhlle_result_dir", 0.13, 1.0, 0.6, 0.0, 0.21, 120.0, 121, 121, 5, -12.0, 12.0, -12.0, 12.0, -2.0, 2.0, 0.6, 25.0, 0.04)
#	generate_hydro_input("vhlle_result_dir", 0.12, 0.85, 0.4, 0.0, 0.21, 145.0, 101, 101, 41, -5.0, 5.0, -5.0, 5.0, -8.0, 8.0, 0.6, 25.0, 0.04)
	print "generate escape config file"
	generate_escape_input("escape_result_dir", 20, 0.21)


	generate_stat_input(   
	m_eta_L = -0.8,
        m_eta_H = 0.8,
        m_pt_L = 0,
        m_pt_H = 100.0,
        dn_eta_L = -5.5,
        dn_eta_H = 5.5,
        dn_eta_bins = 11,
        dn_pt_L = 0.0,
        dn_pt_H = 100.0,
        q_eta_L = -5.6,
        q_eta_H = 5.6,
        q_eta_bins = 7,
        q_pt_L = 0.2,
        q_pt_H = 5.0)


