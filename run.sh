#!/bin/bash
#PBS -q tiny
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=60GB
#PBS -m n
#PBS -e mpiout/
#PBS -o mpiout/
#PBS -N run_mn

export VAR=value
### set the directory of Monte Python code
dir_montepython=
#dir_montepython=${HOME}/work/downloads/montepython_public-2.2.2/

### set the directory of HSC pseudoCl likelihood
likelihoodname=HSC_Y1_pseudoCl_likelihood
dir_likelihood=${dir_montepython}montepython/likelihoods/${likelihoodname}/
if [ ! -d ${dir_likelihood} ]; then
mkdir $dir_likelihood
fi

### working directory
workdir=${dir_montepython}test_code/
if [ ! -d ${workdir} ]; then
mkdir $workdir
fi
cd ${workdir}

HSCdatadir=${workdir}data
if [ ! -d ${HSCdatadir} ]; then
mkdir $HSCdatadir
fi

cp -f __init__.py ${dir_likelihood} 

### number of cores
npro=1

### choose sampler [multinest, MH]
#sampler=multinest
sampler=MH

### multinest parameters
# sampleing efficiency
efr=0.8 
#efr=0.3

#evidence_tolerance
tol=0.5
#tol=0.1

#max iteration number
niter=1000000

#number of live points
nlive=2000

#importance nested sampling
ISflag=True

### the following file is used when running MH sampler
covfile=fiducial_HSC.covmat

### cosmology [free,wmap9,planck15]
#cosmo=free
cosmo=wmap9


#kind: mock,obs
kind=obs

#which spectrum is used: EE, BB, EB
specmode=EE

#cmodel: LCDM,wCDM
cmodel=LCDM

#matter power spectrum: [dmonly, baryon, owls_agn, owls_dm]
matterspec=dmonly

#nuisance parameteers: [fid,noIA,etafix,noshapeerr,nophotozerr,IAonly]
nuisance=fid

#select the method how to derive the photo-z distribution and separate into tomographic bins
pz_method=ephor_ab_rewei_ab
nzmax=299
zcut_method=ephor_ab
if [ $kind = "mock" ]; then
pz_method=mlz
zcut_method=mlz
nzmax=599
fi

#neutrino mass [0,min,vary]
massnu=no

#covariance is fixed or varied: [vary,fix]
covfix=vary

#output covariance in some fixed cosmlogy? [y or n]
# if yes, cosmo set to be wmap9/planck15. Output of the covariance is "cov_fixed.dat", which will be used as the fixed covariance in the following run.
cov_output='n'

#number of tomographic bins [4]
nzbin=4

### write data file for the class ${likelihoodname}
infdata=${dir_likelihood}/${likelihoodname}.data

echo "HSC_Y1_pseudoCl_likelihood.data_directory = '"${HSCdatadir}"'" > $infdata
echo "HSC_Y1_pseudoCl_likelihood.dtype = '"${kind}"'" >> $infdata
echo "HSC_Y1_pseudoCl_likelihood.specmode = '"${specmode}"'" >> $infdata
echo "HSC_Y1_pseudoCl_likelihood.photoz_zcut = '"${zcut_method}"'" >> $infdata
echo "HSC_Y1_pseudoCl_likelihood.photoz_method = '"${pz_method}"'" >> $infdata
echo "HSC_Y1_pseudoCl_likelihood.k_max = 100." >> $infdata
echo "HSC_Y1_pseudoCl_likelihood.nzmax = "${nzmax} >> $infdata
echo "HSC_Y1_pseudoCl_likelihood.nellsmax = 39" >> $infdata
if [ $matterspec = "baryon" ] || [ $matterspec = "owls_agn" ] || [ $matterspec = "owls_dm" ]; then
echo "HSC_Y1_pseudoCl_likelihood.mode = 'hmcode'" >> $infdata
else
echo "HSC_Y1_pseudoCl_likelihood.mode = 'halofit'" >> $infdata
fi
echo "HSC_Y1_pseudoCl_likelihood.zbin_min = [0.30, 0.60, 0.90, 1.20]" >> $infdata
echo "HSC_Y1_pseudoCl_likelihood.zbin_max = [0.60, 0.90, 1.20, 1.50]" >> $infdata
### tomographic bins to be used in the cosmology analysis [1: use, 0:unuse]
echo "HSC_Y1_pseudoCl_likelihood.use_bin = [1, 1, 1, 1]" >> $infdata
### multipole bins to be used in the cosmology analysis [1: use, 0:unuse]
echo "HSC_Y1_pseudoCl_likelihood.bands_to_use = [0, 1, 1, 1, 1, 1, 1, 0, 0]" >> $infdata
case "$covfix" in
    "vary" ) echo "HSC_Y1_pseudoCl_likelihood.cov_recalculate = True" >> $infdata ;;
    "fix" ) echo "HSC_Y1_pseudoCl_likelihood.cov_recalculate = False" >> $infdata ;;
esac
case "$kind" in 
    "mock" ) echo "HSC_Y1_pseudoCl_likelihood.correction_mbias = False" >> $infdata ;;
    "obs" ) echo "HSC_Y1_pseudoCl_likelihood.correction_mbias = True" >> $infdata ;;
esac
if [ $cov_output = 'y' ]; then
echo "HSC_Y1_pseudoCl_likelihood.cov_output = True" >> $infdata
else
echo "HSC_Y1_pseudoCl_likelihood.cov_output = False" >> $infdata
fi

zshiftpara="'zshift_1','zshift_2','zshift_3','zshift_4'"
shapepara="'m_corr', 'a_psf', 'b_psf'"

case "$nuisance" in 
    "fid" ) echo "HSC_Y1_pseudoCl_likelihood.use_nuisance = ["${shapepara}", "${zshiftpara}", 'A_IA', 'exp_IA']" >> $infdata ;;
    "noIA" ) echo "HSC_Y1_pseudoCl_likelihood.use_nuisance = ["${shapepara}", "${zshiftpara}"]" >> $infdata ;;
    "etafix" ) echo "HSC_Y1_pseudoCl_likelihood.use_nuisance = ["${shapepara}", "${zshiftpara}", 'A_IA']" >> $infdata ;;
    "noshapeerr" ) echo "HSC_Y1_pseudoCl_likelihood.use_nuisance = ["${zshiftpara}", 'A_IA', 'exp_IA']" >> $infdata ;;
    "nophotozerr" ) echo "HSC_Y1_pseudoCl_likelihood.use_nuisance = ["${shapepara}", 'A_IA', 'exp_IA']" >> $infdata ;;
    "IAonly" ) echo "HSC_Y1_pseudoCl_likelihood.use_nuisance = ['A_IA', 'exp_IA']" >> $infdata ;;
esac

cat $infdata

### write input parameter file for MCMC samplers
paramfile=input.param
echo "data.experiments=['HSC_Y1_pseudoCl_likelihood']" > $paramfile
if [ $cosmo = "free" ]; then
### default prior
echo "data.parameters['omega_b']      = [0.022,  0.019, 0.026, 0.002,  1, 'cosmo']" >> $paramfile
echo "data.parameters['omega_cdm']    = [0.12,   0.03, 0.7, 0.01,    1, 'cosmo']" >> $paramfile
echo "data.parameters['n_s']          = [0.97,   0.87, 1.07, 0.02,    1, 'cosmo']" >> $paramfile
echo "data.parameters['h']         	= [0.7,     0.6,   0.9,   0.02, 1, 'cosmo']" >> $paramfile
echo "data.parameters['ln10^{10}A_s'] = [3.1,   1.5, 6, 0.1,    1, 'cosmo']" >> $paramfile
elif [ $cosmo = "wmap9" ]; then
echo "data.cosmo_arguments['omega_b'] = 0.02254" >> $paramfile
echo "data.cosmo_arguments['omega_cdm'] = 0.11417" >> $paramfile
echo "data.cosmo_arguments['n_s'] = 0.97" >> $paramfile
echo "data.cosmo_arguments['h'] = 0.6982" >> $paramfile
echo "data.cosmo_arguments['ln10^{10}A_s'] = 3.07" >> $paramfile
elif [ $cosmo = "planck15" ]; then
echo "data.cosmo_arguments['omega_b'] = 0.022256" >> $paramfile
echo "data.cosmo_arguments['omega_cdm'] = 0.11979" >> $paramfile
echo "data.cosmo_arguments['n_s'] = 0.96507" >> $paramfile
echo "data.cosmo_arguments['h'] = 0.6779" >> $paramfile
echo "data.cosmo_arguments['ln10^{10}A_s'] = 3.0891" >> $paramfile
fi

if [ $massnu = "vary" ]; then
echo "data.parameters['m_ncdm']      = [0.06, 0, 1, 0.1, 1, 'cosmo']" >> $paramfile
elif [ $massnu = "min" ]; then
echo "data.cosmo_arguments['m_ncdm'] = 0.06" >> $paramfile
fi

if [ $cmodel = "wCDM" ]; then
echo "data.cosmo_arguments['Omega_Lambda']  = 0" >> $paramfile
echo "data.parameters['w0_fld']    = [-1.0,  -2.0,  -0.3333, 0.1,   1, 'cosmo']" >> $paramfile
fi

if [ $nuisance != "noIA" ] && [ $nuisance != "none" ]; then
echo "data.parameters['A_IA'] = [0., -5., 5., 0.1, 1, 'nuisance']" >> $paramfile
fi
if [ $nuisance != "etafix" ] && [ $nuisance != "noIA" ] && [ $nuisance != "none" ]; then
echo "data.parameters['exp_IA'] = [0., -5., 5., 0.1, 1, 'nuisance']" >> $paramfile
fi
if [ $nuisance != "noshapeerr" ] && [ $nuisance != "none" ] && [ $nuisance != "IAonly" ]; then
echo "data.parameters['m_corr'] = [0.5, 0, 1, 0.1, 1, 'nuisance']" >> $paramfile
echo "data.parameters['a_psf'] = [0.5, 0, 1, 0.1, 1, 'nuisance']" >> $paramfile
echo "data.parameters['b_psf'] = [0.5, 0, 1, 0.1, 1, 'nuisance']" >> $paramfile
fi
if [ $nuisance != "nophotozerr" ] && [ $nuisance != "none" ] && [ $nuisance != "IAonly" ]; then
izbin=0
while [ $izbin -lt $nzbin ];
do
izbin=`expr $izbin + 1`
echo "data.parameters['zshift_"${izbin}"'] = [0.5, 0, 1, 0.1, 1, 'nuisance']" >> $paramfile
done
fi
if [ $matterspec = "baryon" ]; then
echo "data.parameters['c_min'] = [3.13, 1, 5, 0.2, 1, 'cosmo']" >> $paramfile
elif [ $matterspec = "owls_agn" ]; then
echo "data.parameters['c_min'] = [2.32, 2.32, 2.32, 0, 1, 'cosmo']" >> $paramfile
elif [ $matterspec = "owls_dm" ]; then
echo "data.parameters['c_min'] = [3.43, 3.43, 3.43, 0, 1, 'cosmo']" >> $paramfile
echo "data.parameters['eta_0'] = [0.64, 0.64, 0.64, 0, 1, 'cosmo']" >> $paramfile
fi
echo "data.parameters['Omega_m']      = [1, None, None, 0, 1, 'derived']" >> $paramfile
echo "data.parameters['sigma8']       = [1, None, None, 0, 1, 'derived']" >> $paramfile

echo "data.cosmo_arguments['output'] = 'mPk'" >> $paramfile
echo "data.cosmo_arguments['P_k_max_h/Mpc'] = 100." >> $paramfile

echo "data.cosmo_arguments['Omega_k'] = 0." >> $paramfile

if [ $massnu = "min" ] || [ $massnu = "vary" ]; then
echo "data.cosmo_arguments['N_ncdm'] = 1" >> $paramfile
echo "data.cosmo_arguments['deg_ncdm'] = 1" >> $paramfile
echo "data.cosmo_arguments['T_ncdm'] = 0.71611" >> $paramfile
echo "data.cosmo_arguments['N_ur'] = 2.0328" >> $paramfile
fi

echo "data.cosmo_arguments['sBBN file'] = data.path['cosmo']+'/bbn/sBBN.dat'" >> $paramfile
echo "data.cosmo_arguments['k_pivot'] = 0.05" >> $paramfile
echo "data.write_step=1" >> $paramfile

cat $paramfile

if [ ! -d 'chains' ]; then
mkdir chains
fi

dir_chains=chains/HSC_Y1_${kind}

echo $dir_chains
echo $nlive $niter
echo $paramfile

rm -rf ${dir_chains}
mkdir ${dir_chains}

if [ $sampler = 'multinest' ]; then
mpirun -np ${npro} python ${dir_montepython}montepython/MontePython.py run \
    -p ${paramfile} \
    -o ${dir_chains}/ \
    -m NS \
    --NS_importance_nested_sampling ${ISflag} \
    --NS_evidence_tolerance ${tol} \
    --NS_sampling_efficiency ${efr} \
    --NS_n_live_points ${nlive} \
    --NS_max_iter ${niter}
else
python ${dir_montepython}montepython/MontePython.py run \
    -o ${dir_chains} \
    -p ${paramfile} \
    -N 1000000 \
    -c ${covfile} \
    -f 1.5 -j fast
    --update 300
fi

