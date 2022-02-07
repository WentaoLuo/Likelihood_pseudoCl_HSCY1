#!/usr/bin/env python
#
#######################################################################
# likelihood for the HSC WL power spectrum (pseudo-Cl estimator) #
#######################################################################
#
#  ref.: Hikage et al. PASJ 71 (2019) 43, (arXiv:1809.09148)
#
#######################################################################

from montepython.likelihood_class import Likelihood
import os
import sys
import numpy as np
import scipy.interpolate as itp
from scipy.linalg import cholesky, solve_triangular, det
from scipy.stats import norm
from scipy.special import sici
from scipy import integrate
import time
from colossus.cosmology import cosmology
from colossus.lss import mass_function, bias
from colossus.halo import concentration, profile_nfw
from classy import Class

class HSC_Y1_pseudoCl_likelihood(Likelihood):
    
    def __init__(self, path, data, command_line):

        ### initial setup
        Likelihood.__init__(self, path, data, command_line)

        self.need_cosmo_arguments(data, {'output':'mPk'})

        ### read redshift info of tomographic bins
        self.redshift_bins = []
        for index_zbin in xrange(len(self.zbin_min)):
            redshift_bin = '{:.1f}z{:.1f}'.format(self.zbin_min[index_zbin], self.zbin_max[index_zbin])
            self.redshift_bins.append(redshift_bin)

        # number of z-bins
        self.nzbins = len(self.redshift_bins)

        # number of *unique* correlations between z-bins
        self.ntbins = self.nzbins * (self.nzbins + 1) / 2        
        
        ### read which bands to be used for MCMC analysis
        all_bands_to_use = [] 
        bands_to_skip = []

        self.nlbins = len(self.bands_to_use)

        for iband in xrange(self.nlbins):
            bands_to_skip.append(0)

        self.ntbins_use = 0
        for zbin1 in xrange(self.nzbins):
            for zbin2 in xrange(zbin1 + 1): 
                if self.use_bin[zbin1] == 1 and self.use_bin[zbin2] == 1:
                    all_bands_to_use += self.bands_to_use
                    self.ntbins_use += 1
                else:
                    all_bands_to_use += bands_to_skip

        self.indices_for_bands_to_use = np.where(np.asarray(all_bands_to_use) == 1)[0]

        # this is also the number of points in the datavector
        self.ndata = len(self.indices_for_bands_to_use)

        ### read calibration parameters related to PSF error/leakage, multiplicative bias, photo-z
        if self.dtype == "mock":
            self.a_psf_mean = 0.
            self.a_psf_error = 0.
            self.b_psf_mean = 0.
            self.b_psf_error = 0.
        else:
            fname = os.path.join(self.data_directory, 'psf_modelerror.dat')
            self.a_psf_mean, self.a_psf_error = np.loadtxt(fname, unpack=True)
            fname = os.path.join(self.data_directory, 'psf_leakage.dat')
            self.b_psf_mean, self.b_psf_error = np.loadtxt(fname, unpack=True)

        fname = os.path.join(self.data_directory, 'cl_psferr.dat')
        self.Cl_psf_aa, self.Cl_psf_ab, self.Cl_psf_bb = np.loadtxt(fname,unpack=True)

        fname = os.path.join(self.data_directory, '{:}zbin/m_correction.txt'.format(self.nzbins))   
        self.msel_per_zbin = np.loadtxt(fname, usecols=[1])
        self.mzsel_per_zbin = np.loadtxt(fname, usecols=[2])
        self.merror = np.loadtxt(fname, usecols=[3])

        self.zshift_mean_per_zbin = np.zeros((self.nzbins))
        self.zshift_err_per_zbin = np.zeros((self.nzbins))
        fname = os.path.join(self.data_directory, '{:}zbin/zshift_correction.txt'.format(self.nzbins))   
        if self.nzbins == 1:
            self.zshift_mean_per_zbin[0] = np.loadtxt(fname, usecols=[1])
            self.zshift_err_per_zbin[0] = np.loadtxt(fname, usecols=[2])
        else:
            self.zshift_mean_per_zbin = np.loadtxt(fname, usecols=[1])
            self.zshift_err_per_zbin = np.loadtxt(fname, usecols=[2])

        ### read parameters related to 6 disjoint fields (only for the use of cov_NG calculation)
        fname = os.path.join(self.data_directory, 'field_weight.dat')
        self.wei_per_field = np.loadtxt(fname, usecols=[1])
        self.skyarea_per_field = np.loadtxt(fname, usecols=[2])
        self.wei_per_field /= np.sum(self.wei_per_field)

        ### read Cl binning info.
        fname = os.path.join(self.data_directory, 'cl_binning.dat')
        self.ells_wcen, self.ells_min, self.ells_max, self.nmodeeff = np.loadtxt(fname,unpack=True)
        self.nells_bin=len(self.ells_min)

        ### read shot noise spectra directly estimated from randomly-rotating data

        self.band_noise = np.zeros((self.nells_bin,self.nzbins))
        collect_band_noise = []
        for zbin in xrange(self.nzbins):
            fname = os.path.join(self.data_directory, '{:}zbin/{:}/noise_z{:}_{:}.dat'.format(self.nzbins,self.dtype,zbin+1,self.photoz_zcut))
            collect_band_noise = np.loadtxt(fname)
            for iell_bin in xrange(self.nells_bin):
                self.band_noise[iell_bin,zbin] = collect_band_noise[iell_bin]

        # collect BP per zbin and combine into one array
        collect_bp_in_zbins = []
        collect_ells_bp_in_zbins = []
        for zbin1 in xrange(self.nzbins):
            for zbin2 in xrange(zbin1 + 1):
                # zbin2 first in fname!  
                fname = os.path.join(self.data_directory, '{:}zbin/{:}/band_powers_{:}_z{:}xz{:}_{:}.dat'.format(self.nzbins, self.dtype, self.specmode, zbin2 + 1, zbin1 + 1, self.photoz_zcut))
                extracted_band_powers = np.loadtxt(fname)
                collect_bp_in_zbins.append(extracted_band_powers)
                collect_ells_bp_in_zbins.append(self.ells_wcen)
        self.band_powers = np.asarray(collect_bp_in_zbins).flatten()
        self.ells_band_powers = np.asarray(collect_ells_bp_in_zbins).flatten()

        # read noise spectrum covariance directly estimated from randomly-rotating data 
        index_tbin = 0
        self.covnoi_powers = np.zeros((self.nells_bin,self.nells_bin,self.ntbins))
        covnoi_powers_temp = []
        for zbin1 in xrange(self.nzbins):
            for zbin2 in xrange(zbin1 + 1):
                fname = os.path.join(self.data_directory, '{:}zbin/{:}/covnoi_{:}_z{:}xz{:}_{:}.dat'.format(self.nzbins, self.dtype, self.specmode, zbin2 + 1, zbin1 + 1, self.photoz_zcut))
                covnoi_powers_temp = np.loadtxt(fname)
                for iell_bin in xrange(self.nells_bin):
                    for iell_bin2 in xrange(self.nells_bin):
                        self.covnoi_powers[iell_bin,iell_bin2,index_tbin] = covnoi_powers_temp[iell_bin,iell_bin2]
                index_tbin += 1

        # load all the p_z data here and only once per loop-iteration
        self.pz = np.zeros((self.nzmax, self.nzbins), 'float64')
        self.pz_norm = np.zeros(self.nzbins, 'float64')

        for zbin in xrange(self.nzbins):
            redshift_bin = self.redshift_bins[zbin]
            fname = os.path.join(self.data_directory, 'photoz/pz_{:}_{:}.dat'.format(self.photoz_method, redshift_bin))
            print fname
            # z-range is always the same!
            z_hist, n_z_hist = np.loadtxt(fname, unpack=True)

            # I assume that all photo-z have the same z_min, z_max
            self.redshifts = z_hist
            self.pz[:, zbin] = n_z_hist
            # integrate p(z) over z (in view of normalizing it to one)
            dz = self.redshifts[1:] - self.redshifts[:-1]    
            self.pz_norm[zbin] = np.sum(0.5 * (self.pz[1:, zbin] + self.pz[:-1, zbin]) * dz)

        self.z_max = self.redshifts.max()
        for zbin in xrange(self.nzbins):
            self.pz[:, zbin] /= self.pz_norm[zbin]

        self.need_cosmo_arguments(data, {'z_max_pk': self.z_max, 'output': 'mPk', 'non linear': self.mode, 'P_k_max_h/Mpc': self.k_max})
        
        if self.mode == 'hmcode':
            self.cosmo_halofit = Class()
            self.cosmo_owls_dmonly = Class()

        return
      
    def loglkl(self, cosmo, data):

        start_load = time.time()

        self.Omega_m = cosmo.Omega_m() 
        if self.Omega_m > 1:
            return -1e12

        self.Omega_de = cosmo.Omega_Lambda()
        self.Omega_k = 1 - self.Omega_m - self.Omega_de
        self.Omega_b = cosmo.Omega_b()
        self.omega_b = cosmo.omega_b()
        self.omega_cdm = cosmo.Omega0_cdm()*cosmo.h()*cosmo.h()
        self.Omega_nu = cosmo.Omega_nu
        self.sigma8 = cosmo.sigma8()
        self.ns = cosmo.n_s()
        self.small_h = cosmo.h()
        self.H0 = 100*self.small_h 
        c_over_H0 = 2997.9 / self.small_h
        rho_crit = 2.77536627e11 * self.small_h**2
        self.rho_m = rho_crit * self.Omega_m
        self.ln10As=cosmo.get_current_derived_parameters(['ln10^{10}A_s'])['ln10^{10}A_s']

        #photo-z shift parameters:
        zshift_per_zbin = np.zeros(self.nzbins)
        for zbin in xrange(self.nzbins):
            param_name = 'zshift_{:}'.format(zbin + 1)
            if param_name in data.mcmc_parameters:
                rand_z = data.mcmc_parameters[param_name]['current'] * data.mcmc_parameters[param_name]['scale']
                zshift_per_zbin[zbin] = self.zshift_mean_per_zbin[zbin] + self.zshift_err_per_zbin[zbin] * norm.ppf(rand_z)
            else:
                zshift_per_zbin[zbin] = 0.

        #shift photo-z distribution: pz_i (zs) -> pz_i (zs+zshift)
        pz_shifted = np.zeros((self.nzmax, self.nzbins), 'float64')
        for zbin in xrange(self.nzbins):
            interp1d_pz=itp.interp1d(self.redshifts,self.pz[:,zbin])
            for index_z, z in enumerate(self.redshifts):
                old_z = z - zshift_per_zbin[zbin]
                if old_z < self.redshifts.min() or old_z > self.redshifts.max():
                    pz_shifted[index_z,zbin] = 0.
                else:
                    pz_shifted[index_z,zbin]=interp1d_pz(old_z)

        # parameters of PSF modeling error/PSF leakage correction:
        param_name = 'a_psf'
        if param_name in data.mcmc_parameters:
            rand_a = data.mcmc_parameters[param_name]['current'] * data.mcmc_parameters[param_name]['scale']
            a_psf = self.a_psf_mean + self.a_psf_error * norm.ppf(rand_a)
        else:
            a_psf = 0.

        param_name = 'b_psf'
        if param_name in data.mcmc_parameters:
            rand_b = data.mcmc_parameters[param_name]['current'] * data.mcmc_parameters[param_name]['scale']
            b_psf = self.b_psf_mean + self.b_psf_error * norm.ppf(rand_b)
        else:
            b_psf = 0.

        # m-correction:        
        # Errors on m-corrections for different z-bins are correlated, thus one free nuisance "m_corr" is enough
        param_name = 'm_corr'
        if param_name in data.mcmc_parameters:
            rand_m = data.mcmc_parameters[param_name]['current'] * data.mcmc_parameters[param_name]['scale']
            if self.nzbins == 1:
                delta_m_corr = self.merror * norm.ppf(rand_m)
            else:
                delta_m_corr = self.merror[0] * norm.ppf(rand_m)
        else:
            delta_m_corr = 0.

        m_corr_per_zbin = np.zeros(self.nzbins)
        for zbin in xrange(self.nzbins):
            if self.nzbins == 1:
                m_corr_per_zbin[zbin] = self.msel_per_zbin + self.mzsel_per_zbin
            else:
                m_corr_per_zbin[zbin] = self.msel_per_zbin[zbin] + self.mzsel_per_zbin[zbin]

        # needed for IA modelling: 
        if ('A_IA' in data.mcmc_parameters) and ('exp_IA' in data.mcmc_parameters):
            amp_IA = data.mcmc_parameters['A_IA']['current'] * data.mcmc_parameters['A_IA']['scale']
            exp_IA = data.mcmc_parameters['exp_IA']['current'] * data.mcmc_parameters['exp_IA']['scale']
            intrinsic_alignment = True
        elif ('A_IA' in data.mcmc_parameters) and ('exp_IA' not in data.mcmc_parameters):
            amp_IA = data.mcmc_parameters['A_IA']['current'] * data.mcmc_parameters['A_IA']['scale']
            # redshift-scaling is fixed to be 3 if exp_IA is not included:            
            exp_IA = 3.
            intrinsic_alignment = True
        else:
            intrinsic_alignment = False

        # derive the linear growth factor D(z)
        self.linear_growth_rate = np.zeros_like(self.redshifts)

        #print self.redshifts
        for index_z, z in enumerate(self.redshifts):
            try:
                # for CLASS ver >= 2.6:
                self.linear_growth_rate[index_z] = cosmo.scale_independent_growth_factor(z)
            except:
                # my own function from private CLASS modification: 
                self.linear_growth_rate[index_z] = cosmo.growth_factor_at_z(z)
                # normalize to unity at z=0:
            try:
                # for CLASS ver >= 2.6:
                self.linear_growth_rate /= cosmo.scale_independent_growth_factor(0.)
            except:
                # my own function from private CLASS modification:
                self.linear_growth_rate /= cosmo.growth_factor_at_z(0.)

        #> get distances from cosmo-module:
        r, dzdr = cosmo.z_of_r(self.redshifts)

        #> these are the l-nodes for the derivation of the theoretical Cl:
        self.ells = np.logspace(np.log10(self.ells_min[0]), np.log10(self.ells_max[self.nells_bin-1]), self.nellsmax)

        #> After long and extensive testing:
        #> Don't put calls to Class (i.e. cosmo...) into a loop...
        #> before "pk" and the constants were just called at demand below in the code (due to convenience an copy & paste)
        #> which seemed to have been the source for the memory leak...

        if self.mode == 'hmcode':
            current_paramset={'output':'mPk','z_max_pk': self.z_max,'non linear':'halofit','omega_b':self.omega_b,'omega_cdm':self.omega_cdm,'h':self.small_h,'n_s':self.ns,'ln10^{10}A_s':self.ln10As,'P_k_max_h/Mpc': self.k_max}
            try:
                self.cosmo_halofit.set(current_paramset)
                self.cosmo_halofit.compute()
            except:
                print 'HMcode calculation fails:',self.Omega_m,self.Omega_nu,self.small_h,self.ns,self.sigma8
                return -1e12

            current_paramset.update({'non linear':'hmcode','feedback model':'owls_dmonly'})
            #current_paramset.update({'non linear':'hmcode','c_min':'3.13'})
            try:
                self.cosmo_owls_dmonly.set(current_paramset)
                self.cosmo_owls_dmonly.compute()
            except:
                print 'HMcode calculation fails:',self.Omega_m,self.Omega_nu,self.small_h,self.ns,self.sigma8
                return -1e12

        #> Get power spectrum P(k=l/r,z(r)) from cosmological module
        #> this doesn't really have to go into the loop over fields!
        pk = np.zeros((self.nellsmax, self.nzmax),'float64')
        for index_ells in xrange(self.nellsmax):
            for index_z in xrange(self.nzmax):
                k = (self.ells[index_ells] + 0.5) / r[index_z] 
                if k > self.k_max*self.small_h:
                    pk_dm = 0.
                else:
                    pk_dm = cosmo.pk(k, self.redshifts[index_z])

                if self.mode == 'hmcode':
                    pk_dm_model_corr = self.cosmo_halofit.pk(k, self.redshifts[index_z]) / self.cosmo_owls_dmonly.pk(k, self.redshifts[index_z])
                    pk[index_ells, index_z] = pk_dm * pk_dm_model_corr
                    #pk[index_ells, index_z] = pk_dm * pk_dm_model_corr - self.cosmo_halofit.pk(k, self.redshifts[index_z])
                else:
                    pk[index_ells, index_z] = pk_dm

#                if self.baryon_feedback:
#                    if 'A_bary' in data.mcmc_parameters: 
#                        A_bary = data.mcmc_parameters['A_bary']['current'] * data.mcmc_parameters['A_bary']['scale']
#                        pk[index_ells, index_z] = pk_dm * self.baryon_feedback_bias_sqr(k/self.small_h, self.redshifts[index_z], A_bary=A_bary)
#                    else:
#                        pk[index_ells, index_z] = pk_dm * self.baryon_feedback_bias_sqr(k/self.small_h, self.redshifts[index_z])
#                else

        pr = pz_shifted * dzdr[:, np.newaxis]


        #> Compute function g_i(r), that depends on r and z bin
        #> g_i(r) = (3./2.*Omega_m) / (c/H_0)^2 * r(1+z(r)) int_r^+\infty drs P_z(rs) (rs-r)/rs
        self.WL_kernel = np.zeros((self.nzmax, self.nzbins), 'float64')
        for zbin in xrange(self.nzbins):
            for nr in xrange(self.nzmax - 1):
                fun = pr[nr:, zbin] * (r[nr:] - r[nr]) / r[nr:]
                self.WL_kernel[nr, zbin] = np.sum(0.5 * (fun[1:] + fun[:-1]) * (r[nr + 1:] - r[nr:-1]))
                self.WL_kernel[nr, zbin] *= (3. / 2. * self.Omega_m) * r[nr] * (1. + self.redshifts[nr]) / c_over_H0 **2

        #> Start loop over l for computation of C_l^shear
        Cl_GG_integrand = np.zeros((self.nzmax, self.ntbins), 'float64')
        Cl_GG = np.zeros((self.nellsmax, self.ntbins), 'float64')
        Cl = np.zeros((self.nellsmax, self.ntbins), 'float64')

        if intrinsic_alignment:
            Cl_II_integrand = np.zeros_like(Cl_GG_integrand)
            Cl_II = np.zeros_like(Cl_GG)
            Cl_GI_integrand = np.zeros_like(Cl_GG_integrand)
            Cl_GI = np.zeros_like(Cl_GG)

        dr = r[1:] - r[:-1]
        for index_ells in xrange(self.nellsmax):
            #> find Cl_integrand = (g(r) / r)**2 * P(l/r,z(r))
            index_tbin = 0
            for zbin1 in xrange(self.nzbins):
                for zbin2 in xrange(zbin1 + 1):
                    Cl_GG_integrand[:, index_tbin] = self.WL_kernel[:, zbin1] * self.WL_kernel[:, zbin2] / (r[:]**2) * pk[index_ells, :]
                    if intrinsic_alignment:
                        factor_IA = self.get_factor_IA(self.redshifts[:], self.linear_growth_rate[:], amp_IA, exp_IA)
                        Cl_II_integrand[:, index_tbin] = pr[:, zbin1] * pr[:, zbin2] * factor_IA**2 / r[:]**2 * pk[index_ells, :]
                        Cl_GI_integrand[:, index_tbin] = (self.WL_kernel[:, zbin1] * pr[:, zbin2] + self.WL_kernel[:, zbin2] * pr[:, zbin1]) * factor_IA / r[:]**2 * pk[index_ells, :]
                    index_tbin += 1

            #> Integrate over r to get C_l^shear_ij = P_ij(l)
            #> C_l^shear_ii = \sum_0^rmax dr (g_i(r) g_j(r) /r**2) P(k=l/r,z(r))
            index_tbin = 0
            for zbin1 in xrange(self.nzbins):
                for zbin2 in xrange(zbin1 + 1):
                    Cl_GG[index_ells, index_tbin] = np.sum(0.5 * (Cl_GG_integrand[1:, index_tbin] + Cl_GG_integrand[:-1, index_tbin]) * dr)
                    if intrinsic_alignment:
                        Cl_II[index_ells, index_tbin] = np.sum(0.5 * (Cl_II_integrand[1:, index_tbin] + Cl_II_integrand[:-1, index_tbin]) * dr)
                        Cl_GI[index_ells, index_tbin] = np.sum(0.5 * (Cl_GI_integrand[1:, index_tbin] + Cl_GI_integrand[:-1, index_tbin]) * dr)
                        index_tbin += 1

        if self.specmode == 'EE':
            if intrinsic_alignment:
                Cl = Cl_GG + Cl_GI + Cl_II
            else:
                Cl = Cl_GG

        if self.correction_mbias:        
            # add m-correction to theory
            index_tbin = 0
            for zbin1 in xrange(self.nzbins):
                for zbin2 in xrange(zbin1 + 1):
                    m_corr = (1. + m_corr_per_zbin[zbin1]) * (1. + m_corr_per_zbin[zbin2])
                    Cl [:, index_tbin] *= m_corr
                    Cl_GG [:, index_tbin] *= m_corr
                    if intrinsic_alignment:
                        Cl_GI [:, index_tbin] *= m_corr
                        Cl_II [:, index_tbin] *= m_corr
                    index_tbin += 1

        for index_tbin in xrange(self.ntbins):
            Cl [:, index_tbin] *= (1. + delta_m_corr)**2
            Cl_GG [:, index_tbin] *= (1. + delta_m_corr)**2
            if intrinsic_alignment:
                Cl_GI [:, index_tbin] *= (1. + delta_m_corr)**2
                Cl_II [:, index_tbin] *= (1. + delta_m_corr)**2

        # ordering of redshift bins is correct in definition of theory below!
        theory = np.zeros((self.ntbins, self.nells_bin), 'float64')
        theory_GG = np.zeros((self.ntbins, self.nells_bin), 'float64')
        if intrinsic_alignment:
            theory_GI = np.zeros((self.ntbins, self.nells_bin), 'float64')
            theory_II = np.zeros((self.ntbins, self.nells_bin), 'float64')
        cl_psf = np.zeros((self.ntbins, self.nells_bin), 'float64')

        if self.specmode == 'EE':
            nstep = 3
            nells_in_bin = 2**nstep + 1
            index_tbin = 0
            for zbin1 in xrange(self.nzbins):
                for zbin2 in xrange(zbin1 + 1):
                    Cl_sample = Cl[:, index_tbin] 
                    spline_Cl = itp.splrep(self.ells, Cl_sample)
                    Cl_GG_sample = Cl_GG[:, index_tbin] 
                    spline_Cl_GG = itp.splrep(self.ells, Cl_GG_sample)
                    if intrinsic_alignment:
                        Cl_GI_sample = Cl_GI[:, index_tbin] 
                        spline_Cl_GI = itp.splrep(self.ells, Cl_GI_sample)
                        Cl_II_sample = Cl_II[:, index_tbin] 
                        spline_Cl_II = itp.splrep(self.ells, Cl_II_sample)
                    for iell_bin in xrange(self.nells_bin):
                        ells_in_bin = np.linspace(self.ells_min[iell_bin], self.ells_max[iell_bin], nells_in_bin)
                        ells_norm = ells_in_bin * (ells_in_bin + 1) / (2. * np.pi)
                        wei_l = ells_in_bin 
                        D_l = ells_norm * itp.splev(ells_in_bin, spline_Cl)
                        theory[index_tbin,iell_bin] = integrate.simps(D_l*wei_l,ells_in_bin) / integrate.simps(wei_l,ells_in_bin)
                        D_l = ells_norm * itp.splev(ells_in_bin, spline_Cl_GG)
                        theory_GG[index_tbin,iell_bin] = integrate.simps(D_l*wei_l,ells_in_bin) / integrate.simps(wei_l,ells_in_bin)
                        if intrinsic_alignment:
                            D_l = ells_norm * itp.splev(ells_in_bin, spline_Cl_GI)
                            theory_GI[index_tbin,iell_bin] = integrate.simps(D_l*wei_l,ells_in_bin) / integrate.simps(wei_l,ells_in_bin)
                            D_l = ells_norm * itp.splev(ells_in_bin, spline_Cl_II)
                            theory_II[index_tbin,iell_bin] = integrate.simps(D_l*wei_l,ells_in_bin) / integrate.simps(wei_l,ells_in_bin)

                    # PSF modeling error and leakage terms are added to theory                    
                    cl_psf[index_tbin,:] = a_psf * a_psf * self.Cl_psf_aa + 2 * a_psf * b_psf * self.Cl_psf_ab + b_psf * b_psf * self.Cl_psf_bb
                    theory[index_tbin,:] += a_psf * a_psf * self.Cl_psf_aa + 2 * a_psf * b_psf * self.Cl_psf_ab + b_psf * b_psf * self.Cl_psf_bb

                    index_tbin += 1

        band_powers_theory = theory.flatten()
        band_powers_theory_GG = theory_GG.flatten()
        if intrinsic_alignment:
            band_powers_theory_GI = theory_GI.flatten()
            band_powers_theory_II = theory_II.flatten()
        band_powers_cl_psf = cl_psf.flatten()

        if self.cov_recalculate:
            cov = self.cov_GA(Cl)
            #print 'Time for Gaussian covariance:', time.time() - start_load

            if self.specmode == 'EE':
                #!!! low loglike value returns in extremely strange cosmology
                try:
                    cov += self.cov_NG(Cl, cosmo)
                except:
                    print 'cov_NG cannot be calculated due to strange cosmology:',self.Omega_m,self.Omega_nu,self.small_h,self.ns,self.sigma8
                    return -1e12
            # some numpy-magic for slicing:
            cov_sliced = cov[np.ix_(self.indices_for_bands_to_use, self.indices_for_bands_to_use)]
            if self.cov_output:
                corr_sliced = np.zeros((self.ndata,self.ndata))
                for idata1 in xrange(len(cov_sliced)):
                    for idata2 in xrange(len(cov_sliced)):
                        corr_sliced[idata1,idata2] = cov_sliced[idata1,idata2]/np.sqrt(cov_sliced[idata1,idata1]*cov_sliced[idata2,idata2])
                fname = os.path.join(self.data_directory, 'cov_fixed.dat'.format(self.nzbins,self.dtype))
                np.savetxt(fname, cov_sliced.flatten())
                sys.exit()
        else:
            fname = os.path.join(self.data_directory, 'cov_fixed.dat'.format(self.nzbins,self.dtype))
            print('read covariance file:{:}'.format(fname))
            cov_sliced_org = np.loadtxt(fname, usecols=[0])
            cov_sliced_org = np.reshape(cov_sliced_org,(self.nlbins*self.ntbins,self.nlbins*self.ntbins))
            bands_to_use = []
            for zbin1 in xrange(self.nzbins):
                for zbin2 in xrange(zbin1 + 1): 
                    if self.use_bin[zbin1] == 1 and self.use_bin[zbin2] == 1:
                        bands_to_use += self.bands_to_use
                    else:
                        bands_to_use += [0,0,0,0,0,0,0,0,0]
            indices_use = np.where(np.asarray(bands_to_use)==1)[0]
            cov_sliced = cov_sliced_org[np.ix_(indices_use, indices_use)]
        
        difference_vector = self.band_powers - band_powers_theory
        #difference_vector = self.band_powers
        difference_vector = difference_vector[self.indices_for_bands_to_use]
        obs_vector = self.band_powers[self.indices_for_bands_to_use]
        theory_vector = band_powers_theory[self.indices_for_bands_to_use]

        logdet = 2e12
        if np.isinf(band_powers_theory).any() or np.isnan(band_powers_theory).any():
            chi2 = logdet
        else:
            for i in xrange(len(cov_sliced)):
                for j in xrange(len(cov_sliced)):
                    cov_sliced[i,j]=cov_sliced[j,i]
            try:
                cholesky_transform = cholesky(cov_sliced, lower=True)
            except:
                print 'cholesky_transform failed'
                return -1e12
            if self.cov_recalculate:
                try:
                    logdet = 2. * np.log(det(cholesky_transform))
                except:
                    print 'logdet cannot be calculated'
                    return -1e12

            yt = solve_triangular(cholesky_transform, difference_vector, lower=True)
            chi2 = yt.dot(yt)

            yt_obs = solve_triangular(cholesky_transform, obs_vector, lower=True)
            yt_theory = solve_triangular(cholesky_transform, theory_vector, lower=True)

            sn_org = np.sqrt(yt_obs.dot(yt_obs))
            sn_new = yt_obs.dot(yt_theory)/np.sqrt(yt_theory.dot(yt_theory))

        if self.cov_recalculate:
            return -0.5 * chi2 - 0.5 * logdet
        else:
            return -0.5 * chi2


    def cov_GA(self,Cl):
        
        # compute Gaussian term of covariance matrix

        noispec_per_zbin = np.zeros((self.nellsmax,self.nzbins), 'float64')
        for zbin in xrange(self.nzbins):
            spline_band_noise=itp.splrep(self.ells_wcen,self.band_noise[:,zbin])
            for index_ells in xrange(self.nellsmax):
                noispec_per_zbin[index_ells,zbin]=itp.splev(self.ells[index_ells],spline_band_noise)

        cov_ga = np.zeros((self.nellsmax,self.ntbins, self.ntbins), 'float64')
        cov_ga_noi = np.zeros((self.nellsmax,self.ntbins, self.ntbins), 'float64')
        cov_ga_cross = np.zeros((self.nellsmax,self.ntbins, self.ntbins), 'float64')
        for index_ells in xrange(self.nellsmax):
            D_l = self.ells[index_ells] * (self.ells[index_ells] + 1) / (2. * np.pi)
            Cl_tot = np.zeros((self.ntbins), 'float64')
            Cl_sig = np.zeros((self.ntbins), 'float64')
            Cl_noi = np.zeros((self.ntbins), 'float64')
            index_tbin = 0
            for zbin1 in xrange(self.nzbins):
                for zbin2 in xrange(zbin1 + 1):
                    if self.specmode == 'EE':
                        Cl_sig[index_tbin] = Cl[index_ells, index_tbin] * D_l
                        Cl_tot[index_tbin] = Cl_sig[index_tbin]
                    if zbin1 == zbin2:
                        Cl_noi[index_tbin] = noispec_per_zbin[index_ells,zbin1]
                        Cl_tot[index_tbin] += Cl_noi[index_tbin]
                    index_tbin += 1
            index_tbin1 = 0
            for zbin1 in xrange(self.nzbins):
                for zbin2 in xrange(zbin1 + 1):
                    index_tbin2 = 0
                    for zbin3 in xrange(self.nzbins):
                        for zbin4 in xrange(zbin3 + 1):
                            cov_ga[index_ells,index_tbin1,index_tbin2] = Cl_tot[self.zbins2tomo(zbin1,zbin3)] * Cl_tot[self.zbins2tomo(zbin2,zbin4)] + \
                                                            Cl_tot[self.zbins2tomo(zbin1,zbin4)] * Cl_tot[self.zbins2tomo(zbin2,zbin3)] 
                            cov_ga_cross[index_ells,index_tbin1,index_tbin2] = Cl_sig[self.zbins2tomo(zbin1,zbin3)] * Cl_noi[self.zbins2tomo(zbin2,zbin4)] + \
                                                            Cl_noi[self.zbins2tomo(zbin1,zbin3)] * Cl_sig[self.zbins2tomo(zbin2,zbin4)] + \
                                                            Cl_sig[self.zbins2tomo(zbin1,zbin4)] * Cl_noi[self.zbins2tomo(zbin2,zbin3)] + \
                                                            Cl_noi[self.zbins2tomo(zbin1,zbin4)] * Cl_sig[self.zbins2tomo(zbin2,zbin3)] 
                            cov_ga_noi[index_ells,index_tbin1,index_tbin2] = Cl_noi[self.zbins2tomo(zbin1,zbin3)] * Cl_noi[self.zbins2tomo(zbin2,zbin4)] + \
                                                            Cl_noi[self.zbins2tomo(zbin1,zbin4)] * Cl_noi[self.zbins2tomo(zbin2,zbin3)] 
                            index_tbin2 += 1
                    index_tbin1 += 1

        cov_ga_bin = np.zeros((self.nells_bin*self.ntbins, self.nells_bin*self.ntbins), 'float64')
        cov_ga_cross_bin = np.zeros((self.nells_bin*self.ntbins, self.nells_bin*self.ntbins), 'float64')
        cov_ga_noi_bin = np.zeros((self.nells_bin*self.ntbins, self.nells_bin*self.ntbins), 'float64')
        nstep = 3
        nells_in_bin = 2**nstep + 1
        for index_tbin1 in xrange(self.ntbins):
            for index_tbin2 in xrange(self.ntbins):
                cov_sample = cov_ga[:, index_tbin1, index_tbin2] 
                cov_cross_sample = cov_ga_cross[:, index_tbin1, index_tbin2] 
                cov_noi_sample = cov_ga_noi[:, index_tbin1, index_tbin2] 
                spline_cov = itp.splrep(self.ells, cov_sample)
                spline_cov_cross = itp.splrep(self.ells, cov_cross_sample)
                spline_cov_noi = itp.splrep(self.ells, cov_noi_sample)
                for iell_bin in xrange(self.nells_bin):
                    ells_in_bin = np.linspace(self.ells_min[iell_bin], self.ells_max[iell_bin], nells_in_bin)
                    cov_l = itp.splev(ells_in_bin, spline_cov)
                    cov_cross_l = itp.splev(ells_in_bin, spline_cov_cross)
                    cov_noi_l = itp.splev(ells_in_bin, spline_cov_noi)
                    wei_l = ells_in_bin
                    idata1 = iell_bin + index_tbin1 * self.nells_bin
                    idata2 = iell_bin + index_tbin2 * self.nells_bin
                    cov_ga_bin[idata1,idata2] = integrate.simps(cov_l*wei_l,ells_in_bin) / integrate.simps(wei_l,ells_in_bin)
                    cov_ga_bin[idata1,idata2] /= self.nmodeeff[iell_bin]
                    cov_ga_cross_bin[idata1,idata2] = integrate.simps(cov_cross_l*wei_l,ells_in_bin) / integrate.simps(wei_l,ells_in_bin)
                    cov_ga_cross_bin[idata1,idata2] /= self.nmodeeff[iell_bin]
                    cov_ga_noi_bin[idata1,idata2] = integrate.simps(cov_noi_l*wei_l,ells_in_bin) / integrate.simps(wei_l,ells_in_bin)
                    cov_ga_noi_bin[idata1,idata2] /= self.nmodeeff[iell_bin]

        #!!! noise covariance is replaced with the one estimated from randomly-rotating data
        for index_tbin1 in xrange(self.ntbins):
            for index_tbin2 in xrange(self.ntbins):
                for iell_bin in xrange(self.nells_bin):
                    for iell_bin2 in xrange(self.nells_bin):
                        idata1 = iell_bin + index_tbin1 * self.nells_bin
                        idata2 = iell_bin2 + index_tbin2 * self.nells_bin
                        if self.specmode != 'EE' and index_tbin1 == index_tbin2:
                            cov_ga_bin[idata1,idata2] = self.covnoi_powers[iell_bin,iell_bin2,index_tbin1]
                        if self.specmode == 'EE' and index_tbin1 == index_tbin2:
                            cov_ga_bin[idata1,idata2] = cov_ga_bin[idata1,idata2] - cov_ga_noi_bin[idata1,idata2] \
                                                        + self.covnoi_powers[iell_bin,iell_bin2,index_tbin1]

        return cov_ga_bin

    def cov_NG(self, Cl, cosmo):
        
        # compute non-Gaussian terms of covariance matrix
        # Here I use Colossus python code developed by Diemer (2017)
        # halo mass is defined as m200
        # !!! ATTENTION: the unit of length/mass is Mpc/h and Msun/h in this function of cov_NG

        ### assuming flat LCDM csmology
        params_for_colossus={'flat':True,'H0':100.*self.small_h,'Om0':self.Omega_m,'Ob0':self.Omega_b,'sigma8':self.sigma8,'ns':self.ns}
        cosmo_for_colossus=cosmology.setCosmology('myCosmo',params_for_colossus)

        ### for the integral over halo mass
        #logmmin = 8.0
        logmmin = 10.0
        logmmax = 16.0
        #dlogm = 0.1
        dlogm = 0.2
        mass = 10.**np.arange(logmmin,logmmax,dlogm)
        nmmax = len(mass)
        dlnm = np.log(10.**logmmax/10.**logmmin)/nmmax
        mass *= np.exp(dlnm*0.5)

        ### for the integral over redshift 
        zmin=0.01
        zmax=1.5
        #dlnz=0.1
        dlnz=0.2
        zs = np.exp(np.arange(np.log(zmin),np.log(zmax),dlnz) + dlnz * 0.5)

        nzmax = len(zs)
        r, dzdr = cosmo.z_of_r(zs)

        #!!! change unit from Mpc to Mpc/h
        r *= self.small_h
        dzdr /= self.small_h

        # linear matter power spectrum
        logkmin = -4.0
        logkmax = 2.0
        dlogk = 0.1
        k_tab = 10.**np.arange(logkmin,logkmax,dlogk)
        Plin_tab = cosmo_for_colossus.matterPowerSpectrum(k_tab,z=0.)
        logk_tab = np.log(k_tab)
        spline_Plin = itp.splrep(logk_tab, Plin_tab)

        # linear growth rate
        spline_gz = itp.splrep(self.redshifts, self.linear_growth_rate)

        nellsmax = self.nells_bin + 1
        ells = np.logspace(np.log10(self.ells_min[0]), np.log10(self.ells_max[self.nells_bin-1]), nellsmax)
        cov_ng = np.zeros((nellsmax, nellsmax ,self.ntbins, self.ntbins), 'float64')

        fsky = np.zeros(len(self.skyarea_per_field))
        invskyareatot = 0.
        for index_field in xrange(len(self.skyarea_per_field)):
            fsky[index_field] = self.skyarea_per_field[index_field] * np.deg2rad(1.)**2
            invskyareatot += self.wei_per_field[index_field]**2 / fsky[index_field]

        ### table of sigmal2
        s_min = min(1., np.sqrt(fsky.min()) * r[0])
        s_max = max(300., np.sqrt(fsky.max()) * r[nzmax-1])
        ns_tab = 50
        s_tab = np.logspace(np.log10(s_min), np.log10(s_max), ns_tab)
        sigmal2_tab = np.zeros(ns_tab)
        for index_s in xrange(ns_tab):
            sigmal2_tab[index_s] = self.sigmal2_square_field(s_tab[index_s],logk_tab,Plin_tab)
        spline_sigmal2 = itp.splrep(s_tab, sigmal2_tab)

        for index_z in xrange(nzmax):

            lingrowth_z = itp.splev(zs[index_z], spline_gz)

            sigw2 = 0.
            for index_field in xrange(len(self.skyarea_per_field)):
                lsky = np.sqrt(fsky[index_field]) * r[index_z]
                sigw2 += self.wei_per_field[index_field]**2 * itp.splev(lsky, spline_sigmal2)

            sigw2 *= r[index_z]**2 * lingrowth_z**2

            dndlnm = mass_function.massFunction(mass,zs[index_z],mdef='200m',model='tinker08',q_out='dndlnM')
            halobias = bias.haloBias(mass,model='tinker10',z=zs[index_z],mdef='200m')

            c = concentration.concentration(mass,'200m',zs[index_z],model='diemer15')

            rhalo = np.zeros(nmmax)
            for index_m in xrange(nmmax):
                profile = profile_nfw.NFWProfile(M=mass[index_m],mdef='200m',z=zs[index_z],c=c[index_m])
                rhalo[index_m] = profile.RDelta(zs[index_z],'200m') * 1e-3   ## change unit from kpc/h to Mpc/h

            wk = np.zeros(self.nzbins)
            for zbin in xrange(self.nzbins):
                interp1d_WL_kernel = itp.interp1d(self.redshifts,self.WL_kernel[:,zbin])
                if zs[index_z] < self.redshifts[0]:
                    wk[zbin] = 0.
                else:
                    wk[zbin] = interp1d_WL_kernel(zs[index_z]) / self.small_h

            wkpro = np.zeros(self.ntbins)
            index_tbin = 0
            for zbin1 in xrange(self.nzbins):
                for zbin2 in xrange(zbin1+1):
                    wkpro[index_tbin] = wk[zbin1] * wk[zbin2]
                    index_tbin += 1

            common_fac = zs[index_z] * dlnz / dzdr[index_z] / r[index_z]**6

            for index_ells1 in xrange(nellsmax):
                D1_l = ells[index_ells1] * (ells[index_ells1] + 1) / (2. * np.pi)
                for index_ells2 in xrange(index_ells1+1):
                    D2_l = ells[index_ells2] * (ells[index_ells2] + 1) / (2. * np.pi)

                    k1 = (ells[index_ells1] + 0.5) / r[index_z] 
                    k2 = (ells[index_ells2] + 0.5) / r[index_z] 

                    pklz = np.zeros(2)
                    pklz[0]=itp.splev(np.log(k1), spline_Plin) * lingrowth_z**2
                    pklz[1]=itp.splev(np.log(k2), spline_Plin) * lingrowth_z**2

                    km = np.zeros((2))
                    clmm = np.zeros((2))
                    clm = np.zeros((2))
                    cl1h = 0.

                    kmout = np.zeros(nmmax)
                    for index_m in xrange(nmmax):
                        km[0] = self.winnfw(k1,rhalo[index_m],c[index_m])
                        km[1] = self.winnfw(k2,rhalo[index_m],c[index_m])
                        km *= mass[index_m] / (self.rho_m / self.small_h**2 )
                        kmout[index_m] = km[0]
                        clmm += km * km * halobias[index_m] * dndlnm[index_m] * dlnm
                        clm += km * halobias[index_m] * dndlnm[index_m] * dlnm
                        cl1h += (km[0] * km[1])**2 * dndlnm[index_m] * dlnm

                    for index_tbin1 in xrange(self.ntbins):
                        for index_tbin2 in xrange(self.ntbins):
                            wkpro2 = wkpro[index_tbin1] * wkpro[index_tbin2]
                            cov_HSV = clmm[0] * clmm[1] * wkpro2 * sigw2 * common_fac
                            cov_BC = (68. / 21.)**2 * (clm[0] * clm[1])**2 * pklz[0] * pklz[1] * wkpro2 * sigw2 * common_fac
                            cov_HSVBC = (68. / 21.) * (clm[0]**2 * clmm[1] * pklz[0] + clm[1]**2 * clmm[0] * pklz[1]) * wkpro2 * sigw2 * common_fac
                            cov_cNG = cl1h * wkpro2 * common_fac * invskyareatot
                            cov_tot = cov_HSV + cov_BC + cov_HSVBC + cov_cNG
                            cov_ng[index_ells1,index_ells2,index_tbin1,index_tbin2] += cov_tot * D1_l * D2_l
                        
        for index_tbin1 in xrange(self.ntbins):
            for index_tbin2 in xrange(self.ntbins):
                for index_ells2 in xrange(nellsmax):
                    for index_ells1 in xrange(index_ells2):
                        cov_ng[index_ells1,index_ells2,index_tbin1,index_tbin2] = cov_ng[index_ells2,index_ells1,index_tbin1,index_tbin2] 

        cov_ng_bin = np.zeros((self.nells_bin*self.ntbins, self.nells_bin*self.ntbins), 'float64') 

        nstep = 2
        nells_in_bin = 2**nstep + 1
        for index_tbin1 in xrange(self.ntbins):
            for index_tbin2 in xrange(self.ntbins):
                cov_sample = cov_ng[:, :, index_tbin1, index_tbin2] 
                interp2d_cov = itp.interp2d(ells, ells, cov_sample, kind='cubic')
                for iell_bin1 in xrange(self.nells_bin):
                    ells_in_bin1 = np.linspace(self.ells_min[iell_bin1], self.ells_max[iell_bin1], nells_in_bin)
                    ells_in_bin_mid1 = 0.5 * (ells_in_bin1[1:] + ells_in_bin1[:-1])
                    wei_l1 = ells_in_bin_mid1 / np.sum(ells_in_bin_mid1)
                    for iell_bin2 in xrange(self.nells_bin):
                        ells_in_bin2 = np.linspace(self.ells_min[iell_bin2], self.ells_max[iell_bin2], nells_in_bin)
                        ells_in_bin_mid2 = 0.5 * (ells_in_bin2[1:] + ells_in_bin2[:-1])
                        wei_l2 = ells_in_bin_mid2 / np.sum(ells_in_bin_mid2)
                        wbin = wei_l1 * wei_l2[:,np.newaxis]
                        cov_int = interp2d_cov(ells_in_bin_mid1,ells_in_bin_mid2)
                        idata1 = iell_bin1 + index_tbin1 * self.nells_bin
                        idata2 = iell_bin2 + index_tbin2 * self.nells_bin
                        cov_ng_bin[idata1,idata2] = np.sum(cov_int*wbin)

        return cov_ng_bin

    def sigmal2_square_field(self,s,logk_tab,Plin_tab):

        sigmal2 = 0.
        num_k = len(logk_tab)
        dlogk = logk_tab[1:]-logk_tab[:-1]
        k = np.exp(logk_tab)
        p = k * s / (2. * np.pi)
        dsigdlnk = k**2 / (2. * np.pi) * Plin_tab * self.win_sq(p)**2
        sigmal2 = np.sum(0.5*(dsigdlnk[1:]+dsigdlnk[:-1])*dlogk)

        return sigmal2

    def sigmal2_rec_field(self,sx,sy,logk_tab,Plin_tab):

        sigmal2 = 0.
        num_k = len(logk_tab)
        dlogk = logk_tab[1:]-logk_tab[:-1]
        k = np.exp(logk_tab)
        px = k * sx / (2. * np.pi)
        py = k * sy / (2. * np.pi)
        dsigdlnk = k**2 / (2. * np.pi) * Plin_tab * self.win_rec(px,py)**2
        sigmal2 = np.sum(0.5*(dsigdlnk[1:]+dsigdlnk[:-1])*dlogk)

        return sigmal2

    def winnfw(self,k,rv,c):
        
        ks = k * rv / c

        si1,ci1 = sici(ks)
        si2,ci2 = sici((1.+c)*ks)

        p1 = np.cos(ks) * (ci2 - ci1)
        p2 = np.sin(ks) * (si2 - si1)
        p3 = np.sin(ks*c) / (ks * (1. + c))

        norm = np.log(1. + c) - c / (1. + c)

        return (p1+p2-p3)/norm

    def win_sq(self,p):

        tmin=0.
        tmax = 0.5 * np.pi
        imax = 50
        dt = (tmax - tmin) / imax

        win = 0.
        for i in xrange(imax):
            t = tmin + dt * (i + 0.5)
            x = p * np.cos(t)
            y = p * np.sin(t)
            win += np.sinc(x)*np.sinc(y) / imax

        return win

    def win_rec(self,px,py):

        imax = 50
        pmax = 0.5 * np.pi
        dp = pmax / imax

        win = 0.
        for i in xrange(imax):
            p = dp * (i + 0.5)
            x = px * np.cos(p)
            y = py * np.sin(p)
            win += np.sinc(x) * np.sinc(y) * dp / pmax

        return win

    def zbins2tomo(self,izbin1,izbin2):

        if izbin1 >= izbin2:
            return (izbin1 * (izbin1 + 1)) / 2 + izbin2
        else:
            return (izbin2 * (izbin2 + 1)) / 2 + izbin1

    
    def baryon_feedback_bias_sqr(self, k, z, A_bary=1.):
        """
        
        Fitting formula for baryon feedback following equation 10 and Table 2 from J. Harnois-Deraps et al. 2014 (arXiv.1407.4301)
        
        """

        # k is expected in h/Mpc and is divided in log by this unit...
        x = np.log10(k)
        
        a = 1. / (1. + z)
        a_sqr = a * a
        
        constant = {'AGN':   {'A2': -0.11900, 'B2':  0.1300, 'C2':  0.6000, 'D2':  0.002110, 'E2': -2.0600, 
                              'A1':  0.30800, 'B1': -0.6600, 'C1': -0.7600, 'D1': -0.002950, 'E1':  1.8400,
                              'A0':  0.15000, 'B0':  1.2200, 'C0':  1.3800, 'D0':  0.001300, 'E0':  3.5700},
                    'REF':   {'A2': -0.05880, 'B2': -0.2510, 'C2': -0.9340, 'D2': -0.004540, 'E2':  0.8580, 
                              'A1':  0.07280, 'B1':  0.0381, 'C1':  1.0600, 'D1':  0.006520, 'E1': -1.7900, 
                              'A0':  0.00972, 'B0':  1.1200, 'C0':  0.7500, 'D0': -0.000196, 'E0':  4.5400},
                    'DBLIM': {'A2': -0.29500, 'B2': -0.9890, 'C2': -0.0143, 'D2':  0.001990, 'E2': -0.8250,
                              'A1':  0.49000, 'B1':  0.6420, 'C1': -0.0594, 'D1': -0.002350, 'E1': -0.0611,
                              'A0': -0.01660, 'B0':  1.0500, 'C0':  1.3000, 'D0':  0.001200, 'E0':  4.4800}}
        
        A_z = constant[self.baryon_model]['A2']*a_sqr+constant[self.baryon_model]['A1']*a+constant[self.baryon_model]['A0']
        B_z = constant[self.baryon_model]['B2']*a_sqr+constant[self.baryon_model]['B1']*a+constant[self.baryon_model]['B0']
        C_z = constant[self.baryon_model]['C2']*a_sqr+constant[self.baryon_model]['C1']*a+constant[self.baryon_model]['C0']
        D_z = constant[self.baryon_model]['D2']*a_sqr+constant[self.baryon_model]['D1']*a+constant[self.baryon_model]['D0']
        E_z = constant[self.baryon_model]['E2']*a_sqr+constant[self.baryon_model]['E1']*a+constant[self.baryon_model]['E0']
        
        bias_sqr = 1. - A_bary * (A_z * np.exp((B_z * x - C_z)**3) - D_z * x * np.exp(E_z * x))
        
        return bias_sqr
    
    def get_factor_IA(self, z, linear_growth_rate, amplitude, exponent):

        const = 5e-14 / self.small_h**2 # in Mpc^3 / M_sol
        z0 = 0.62

        factor = -1. * amplitude * const * self.rho_m / linear_growth_rate * ((1. + z) / (1. + z0))**exponent

        return factor         
    
    
    
    
