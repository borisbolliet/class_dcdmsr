
import os
# os.environ['OPENBLAS_NUM_THREADS'] = '1'
#os.system("export OMP_NUM_THREADS=8")
import argparse
import numpy as np
import subprocess
import multiprocessing
import re
import functools
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from datetime import datetime
from datetime import date
import sys
import yaml
from getdist import loadMCSamples, MCSamples
from getdist import plots

from uncertainties import ufloat
from uncertainties.umath import *


path_to_sz_projects = "/Users/boris/Work/CLASS-SZ/SO-SZ/"
path_to_anaconda3 = '/Users/boris/anaconda3/'

chain_dir = '/Users/boris/Work/CLASS-SZ/SO-SZ/class_T0/mcmcs/mcmc_chains/'


path_to_figures = '/Users/boris/Work/CLASS-SZ/SO-SZ/class_T0/mcmcs/'
# path_to_getdist = path_to_anaconda3 + 'bin/getdist'
path_to_getdist = 'getdist'
###############################################


path_to_chains = []
final_chain_dir_list = []


# here put the  full path to chain excluding 'txt' extension
final_chain_dir_list.append(chain_dir+'/CLASS2p8_planck2018_DCDMSRNeff_logpriors_lowTEB_plikHM_TTTEEE_SH0ES20_FIRAS_BAO_DR7_DR12_6dF/CLASS2p8_planck2018_DCDMSRNeff_logpriors_lowTEB_plikHM_TTTEEE_SH0ES20_FIRAS_BAO_DR7_DR12_6dF')
# here put the  full path to chain directory
path_to_chains.append(chain_dir+'/CLASS2p8_planck2018_DCDMSRNeff_logpriors_lowTEB_plikHM_TTTEEE_SH0ES20_FIRAS_BAO_DR7_DR12_6dF/')


# here put the  full path to chain excluding 'txt' extension
final_chain_dir_list.append(chain_dir+'/CLASS2p8_planck2018_DCDMSR_logpriors_lowTEB_plikHM_TTTEEE_FIRAS/CLASS2p8_planck2018_DCDMSR_logpriors_lowTEB_plikHM_TTTEEE_FIRAS')
# here put the  full path to chain directory
path_to_chains.append(chain_dir+'/CLASS2p8_planck2018_DCDMSR_logpriors_lowTEB_plikHM_TTTEEE_FIRAS/')


# here put the  full path to chain excluding 'txt' extension
# /Users/boris/Work/CLASS-SZ/SO-SZ/class_T0/mcmcs/mcmc_chains/CLASS2p8_planck2018_DCDMSR_logpriors_lowTEB_plikHM_TTTEEE_SH0ES20_FIRAS/CLASS2p8_planck2018_DCDMSR_logpriors_lowTEB_plikHM_TTTEEE_SH0ES20_FIRAS.updated.yaml
final_chain_dir_list.append(chain_dir+'/CLASS2p8_planck2018_DCDMSR_logpriors_lowTEB_plikHM_TTTEEE_SH0ES20_FIRAS/CLASS2p8_planck2018_DCDMSR_logpriors_lowTEB_plikHM_TTTEEE_SH0ES20_FIRAS')
# here put the  full path to chain directory
path_to_chains.append(chain_dir+'/CLASS2p8_planck2018_DCDMSR_logpriors_lowTEB_plikHM_TTTEEE_SH0ES20_FIRAS/')






def run(args):
    for i in range(len(final_chain_dir_list)):
        #i = 0
        print('plotting best-fitting results')
        #step_previous_run = step
        #path_to_previous_run = path_to_chains + "/sz_ps_completeness_analysis_" + final_chain_dir_list[i]
        #path_to_previous_run = path_to_chains + "/sz_ps_completeness_analysis_" + final_chain_dir_list[i]
        #path_to_previous_run = path_to_chains + "/sz_ps_completeness_analysis_" + final_chain_dir_list[i]
        chains_previous_run =  final_chain_dir_list[i]
        chain_name = final_chain_dir_list[i]

        os.chdir(path_to_chains[i])
        print('running getdist on chains ' + final_chain_dir_list[i])
        str_cmd_subprocess = ["nice","-n","19",path_to_getdist,final_chain_dir_list[i],"--ignore_rows","0.1"]
        #UNCOMMENT!

        subprocess.call(str_cmd_subprocess)
        # exit(0)
        path_to_likestats = final_chain_dir_list[i] + '.likestats'
        print(path_to_likestats)
        #path_to_covmat = chains_previous_run + '.covmat'
        with open(path_to_likestats) as f:
            for line in f:
                x = line.strip()
                if x:
                    l = re.split('\s',x)
                    l = [e for e in l if e]
                    if l[0]=="H0":
                        H0 = float(l[1])
                        print("best-fitting: H0 = %.4e"%H0)
                    if l[0]=="n_s":
                        n_s = float(l[1])
                        print("best-fitting: n_s = %.4e"%n_s)
                    if l[0]=="omega_b":
                        omega_b = float(l[1])
                        print("best-fitting: omega_b = %.4e"%omega_b)
                    if l[0]=="omega_cdm":
                        omega_cdm = float(l[1])
                        print("best-fitting: omega_cdm = %.4e"%omega_cdm)
                    if l[0]=="chi2*":
                        chi2 = float(l[1])
                        print("best-fitting: chi2*= %.4e"%chi2)
        path_to_margestats = chains_previous_run + '.margestats'
        flag_params = 0
        path_to_yaml = chains_previous_run + '.updated.yaml'
        with open(path_to_yaml, 'r') as stream:
            try:
                parsed_yaml=yaml.safe_load(stream)
                #print(parsed_yaml)
            except yaml.YAMLError as exc:
                print(exc)
        paramlist = parsed_yaml['theory']['classy']['input_params']+parsed_yaml['theory']['classy']['output_params']
        str_first_param = paramlist[0]
        #create a dictionnary to store mean and stddev values:
        marg_results = {}
        with open(path_to_margestats) as f:
            for line in f:
                x = line.strip()
                if x:
                    l = re.split('\s',x)
                    l = [e for e in l if e]
                    if l[0]== str_first_param:
                        flag_params = 1
                    if flag_params == 1:
                        mean = float(l[1])
                        stddev = float(l[2])
                        print("mean +/- stddev: %s"%l[0] + " = %.4f +/- %.4f"%(mean,stddev))
                        marg_results[l[0]] = {}
                        marg_results[l[0]]['mean'] = mean
                        marg_results[l[0]]['stddev'] = stddev

    #exit(0)

    print('running getdist on chains ')

    chains = []
    for i in range(len(final_chain_dir_list)):
        #path_to_previous_run = path_to_chains + "/sz_ps_completeness_analysis_" + final_chain_dir_list[i]
        chains.append(final_chain_dir_list[i])
    # chains = [
    # chains_previous_run
    #     ]

    all_samples = []
    T_cmb_firas = 2.7255

    for i in range (0,len(chains)):
        print(chains[i])
        readsamps = loadMCSamples(chains[i],settings={'ignore_rows':0.1})
        p = readsamps.getParams()
        try:
            readsamps.addDerived(p.H0*np.power(p.T_cmb/T_cmb_firas,1.)
                                 ,name = 'H0_hat',
                          label = r'\hat{H}_0')
            readsamps.addDerived(p.Omega_m*np.power(p.T_cmb/T_cmb_firas,-4.5)
                                 ,name = 'Omega_m_hat',
                          label = r'\hat{\Omega}_\mathrm{m}')
            readsamps.addDerived(p.sigma8*np.power(p.T_cmb/T_cmb_firas,1.2)
                                 ,name = 'sigma8_hat',
                          label = r'\hat{\sigma}_8')
        except:
            print('couldnt add param')

        #     readsamps.addDerived((p.sigma8/0.8)*(((p.Omega_m))/0.3)**(0.35)*(((p.B))/1.25)**(-0.35)*(p.H0/70.)**(-0.2)
        #                          ,name = 'F',
        #                   label = r'(\frac{\sigma_8}{0.8})(\frac{\Omega_m}{0.3})^{^{0.35}}(\frac{B}{1.25})^{^{-0.35}}h_{70}^{^{-0.2}}')

        readsamps.ranges.setRange('log10_Gamma_dcdm',[2., 7])
        samples = readsamps
        samples.updateBaseStatistics()
        all_samples.append(samples)

        Fs = readsamps.getInlineLatex('H0_hat',limit=1)
        print("H0_hat : ",Fs)
        Fs = readsamps.getInlineLatex('Omega_m_hat',limit=1)
        print("Omega_m_hat : ",Fs)
        Fs = readsamps.getInlineLatex('sigma8_hat',limit=1)
        print("sigma8_hat : ",Fs)
        # Fs = Fs.split('=')
        # Fs = Fs[1].split('\pm')
        # Fs= ufloat(float(Fs[0]), float(Fs[1]))
        # print(final_chain_dir_list[i])
        # print("H0_hat = {:.3u}".format(Fs))

    g = plots.getSubplotPlotter()
    g.settings.fig_width_inch = 12

    g.settings.axes_fontsize = 10
    g.settings.lab_fontsize =13
    g.settings.legend_fontsize = 14
    g.settings.alpha_filled_add=0.9
    g.settings.colorbar_label_pad = 20.
    g.settings.figure_legend_frame = False

    sample_list = []
    for s in range(len(all_samples)):
        sample_list.append(all_samples[s])
    g.triangle_plot(sample_list,
                    [ 'Neff','log10_Gamma_dcdm','log10_omega_ini_dcdm_hat','T_cmb',
                    'omega_b_hat',#'omega_b',
                    'omega_cdm_hat',#'omega_cdm',
                    'logA_hat','n_s','theta_s_1e2','H0','H0_hat','sigma8','sigma8_hat','Omega_m','Omega_m_hat','omegamh2','s8omegamp5_norm'],
                    filled=True,
                    legend_labels=[ 'CMB+FIRAS+SHOES+Neff','CMB+FIRAS', 'CMB+FIRAS+SHOES'],
                    legend_loc='upper right',
                    colors = ['blue','red','green','red'],
                    line_args=[{'lw':'1','color':'blue'},{'lw':'1','color':'red'},{'lw':'1','color':'green'},{'lw':'1','color':'r'}]
                   )


    g.export(path_to_figures+'contours.pdf')



def main():
    parser=argparse.ArgumentParser(description="run mcmc")
    #parser.add_argument("-run_mode",help="mcmc or plot or contours",dest="run_mode",type=str,required=True)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
	main()
