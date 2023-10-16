## --------------
## use: python plot_and_compute_fractions_checkCuts.py processo operatore valore_operatore
## -------------

import argparse

import os

import sys
import numpy as np
from hist import Hist
from lhereader import LHEReader
from matplotlib import pyplot as plt
import json
from cycler import cycler
import gzip
import shutil

def plot(histograms, processo, oppe, valu, outdir):
    '''Plots all histograms. No need to change.'''

    for observable, histogram in histograms.items():
        plt.gcf().clf()
        histogram.plot()
        plt.gca().axvline(x=1100., color='Red')  # line of the mass cut
        plt.gcf().text(0.65, 0.9, oppe + ' ' + valu)  # print on plot the operator and the value
        plt.gcf().savefig(os.path.join(outdir, 'plot_' + processo + '_' + oppe + '_' + valu +'.pdf'))
        plt.gcf().savefig(os.path.join(outdir, 'plot_' + processo + '_' + oppe + '_' + valu +'.png'))

        

def setup_histograms():
    '''Histogram initialization. Add new histos here.'''
    
    # Bin edges for each observable
    # TODO: Add your new observables and binnings here
    bins ={
#        'mass_zzh' : np.linspace(1000,5000,50), # final state zzh
        'mass_wzh' : np.linspace(1000,5000,50), # final state w+zh
#        'mass_zhh' : np.linspace(1000,5000,50), # final state zhh
    } 

    # No need to change this part
    histograms = { 
                    observable : (
                                    Hist.new
                                    .Var(binning, name=observable, label=observable)
                                    .Int64()
                                )
                    for observable, binning in bins.items()
    }

    return histograms

def analyze(processo, oppe, valu, outdir):
    '''Event loop + histogram filling'''

    ### --- gg->zzh
#    lhe_file = os.path.join('/afs', 'cern.ch', 'user', 'a', 'acappati', 'work', 'ZZH', '220219_process1', oppe, 'MG5_aMC_v2_7_3_py3', processo, 'Events', 'run_' + oppe + '_' + valu + '_cuts', 'unweighted_events.lhe') # process 1

    ### --- pp->zhh
#    lhe_file = os.path.join('/afs', 'cern.ch', 'user', 'a', 'acappati', 'work', 'ZZH', '220219_process3', oppe, 'MG5_aMC_v2_7_3_py3', processo, 'Events', 'run_' + oppe + '_' + valu + '_cuts', 'unweighted_events.lhe') # process 3

    ### --- pp->zzh
#    lhe_file = os.path.join('/afs', 'cern.ch', 'user', 'a', 'acappati', 'work', 'ZZH', '231010_v2', 'MG5_aMC_v2_7_3_py3', processo, 'Events', 'run_' + oppe + '_' + valu + '_cuts', 'unweighted_events.lhe') 

    ### --- pp->w+zh
    lhe_file = os.path.join('/afs', 'cern.ch', 'user', 'a', 'acappati', 'work', 'ZZH', '231010_v2', 'MG5_aMC_v2_7_3_py3', processo, 'Events', 'run_' + oppe + '_' + valu + '_cuts', 'unweighted_events.lhe') 

    lhe_file_gz = lhe_file + '.gz'

    # check if gzipped file exists
    if os.path.isfile(lhe_file_gz):
        # open the gzip (read+byte mode) as file_in context
        # open the recipient unzipped file (write+byte mode) as file_out context
        with gzip.open(lhe_file_gz, 'rb') as file_in, open(lhe_file, 'wb') as file_out:
            # copy content of zipped file to recipient
            shutil.copyfileobj(file_in, file_out)

    # check if unzipped file exists
    if os.path.isfile(lhe_file):
        reader = LHEReader(lhe_file)
    else:
        # throw error
        raise FileNotFoundError(f'{lhe_file} not found!')

    histograms = setup_histograms()
    # --- dictionary for results
    limit_list = {
       1000000. : 0 , 
       1200. : 0 ,
       1400. : 0 ,
       1600. : 0 ,
       1800. : 0 ,
       2000. : 0 ,
       2500. : 0 ,
       3000. : 0 ,
       3500. : 0 ,
       4000. : 0 ,
       5000. : 0
    }

    for event in reader:
        # Find tops
        tops = filter(
            lambda x: abs(x.pdgid) in (23,24,25),
            event.particles
        )
        jets = filter(
            lambda x1: abs(x1.pdgid) in (1,2,3,4,5) and x1.status > 0,
            event.particles
        )

        # Sum over top four-momenta in the event
        combined_p4 = None
        for p4 in map(lambda x: x.p4(), tops):
            if combined_p4:
                combined_p4 += p4
            else:
                combined_p4 = p4

        for i_limit in limit_list.keys():
            if combined_p4.mass < i_limit: 
                limit_list[i_limit] += 1
     
        combined_p42 = None
        for p42 in map(lambda x1: x1.p4(), jets):
            if combined_p42:
                combined_p42 += p42
            else:
                combined_p42 = p42

        # --- save file with fit results
        outfile = os.path.join(outdir, 'fractions_' + processo + '_' + oppe + '_' + valu + '.json')
        with open(outfile,'w') as f:
                json.dump(limit_list,f)

        # mass histogram
#        histograms['mass_zzh'].fill(combined_p4.mass, weight=1.) # process 1
        histograms['mass_wzh'].fill(combined_p4.mass, weight=1.)
#        histograms['mass_zhh'].fill(combined_p4.mass, weight=1.) # process 3

    return histograms

def argparser(description):

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('process', help='MADGraph output directory')
    parser.add_argument('operator', help='the name of operator')
    parser.add_argument('value', help='the value of the operator')

    args = parser.parse_args()

    return args

def main():

    description = \
    """
    script to check cuts of the madgraph process
    """
    args = argparser(description)

    processo = args.process
    oppe = args.operator
    valu = args.value


    outdir = './plotsAndFractions_231010/'
    os.makedirs(outdir, exist_ok=True)

    sys.stderr.write('Opening file and forming histo...')
    histograms = analyze(processo, oppe, valu, outdir)
    sys.stderr.write('\nHistogram has been filled.')
    sys.stderr.write('\nSaving plot...')
    plot(histograms, processo, oppe, valu, outdir)
    sys.stderr.write(f'\nPlot(s) saved in {outdir}.\n')

if __name__=="__main__":
    main()
    sys.exit(0)
