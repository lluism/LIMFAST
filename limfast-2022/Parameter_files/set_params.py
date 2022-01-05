# Richard H. Mebane
# set_params.py
# February 5, 2019
# Python class for manipulating 21cmFAST parameter files

import numpy as np

# PARAM DESCRIPTIONS

# if a file name is empty, then use power law params
# if 


def set_defaults():
	pf = \
	{
	"USER_SOURCE_DEF": 0, # if this is 1, must write own C functions for sources as prototyped
	"N_GAMMA_UV": 5000,
	"N_GAMMA_UV_FILE": "-1",
	"STELLAR_BARYON_FRAC": 0.05,
	"STELLAR_BARYON_PL": 0.5,
	"STELLAR_BARYON_FRAC_FILE": "-1",
	"ESC_FRAC": 0.1,
	"ESC_PL": -0.5,
	"ESC_FRAC_FILE": "-1",
	"M_TURNOVER": 5.0e8,
	"USE_GENERAL_SOURCES": 0,
	}
	return pf

# given a parameter dict, write it to file
# fname doesn't include extenstion!
def write_pf(pf, fname):
	data = []
	data.append(pf['USER_SOURCE_DEF'])
	data.append(pf['N_GAMMA_UV'])
	data.append(pf['STELLAR_BARYON_FRAC'])
	data.append(pf['STELLAR_BARYON_PL'])
	data.append(pf['ESC_FRAC'])
	data.append(pf['ESC_PL'])
	data.append(pf['M_TURNOVER'])
	data.append(pf['N_GAMMA_UV_FILE'])
	data.append(pf['STELLAR_BARYON_FRAC_FILE'])
	data.append(pf['ESC_FRAC_FILE'])
	data.append(pf['USE_GENERAL_SOURCES'])
	np.savetxt(fname + ".dat", np.c_[data], fmt="%s")



	

pf = set_defaults()
write_pf(pf, "test")
