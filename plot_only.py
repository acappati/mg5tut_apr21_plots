# run with 
# python3 plot_only.py processo operatore valore_operatore_1(SM) valore_operatore_2 xsec1_inPb xsec2_inPb label 

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

def plot(data1,data2,oppe,valu2,xs,xs2,label):
#def plot(data1, data2,oppe,valu,label):
    '''Plots all histograms. No need to change.'''

    outdir = './plots_vsSM_220414_200kevt/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    hist_name = label + '_mass'
    xaxis_name = "$m_{" + label + "}$ [GeV]"

    nbins = 30

    fig = plt.figure(figsize=(8.,6.))

    # print(data1)
    # print(data2)

    weight1 = np.full_like(data1, (float(xs)*1000)/float(nbins)) # weight: xsec [pb] transformed in fb divided per number of bins
    weight2 = np.full_like(data2, (float(xs2)*1000)/float(nbins))

    print(weight1)
    print(weight2)
    print(xs)
    print(xs2)

    # Plot two distributions on the same plot
    ax1 = fig.add_subplot(3,1,(1,2)) # nrows, ncols, index
    ax1.set_ylabel("cross section [fb/100 GeV]")
    ax1.set_xticks(np.arange(550,3550+300,300)) #set ticks

    val_of_bins_data2, edges_of_bins_data2, patches_data2 = plt.hist(data2, nbins, range=(550,3550), weights=weight2, histtype='step', label="$f_{" + oppe[1:3] + "}/\Lambda^4 =" + valu2 + "$ TeV$^{-4}$")
    val_of_bins_data1, edges_of_bins_data1, patches_data1 = plt.hist(data1, nbins, range=(550,3550), weights=weight1, histtype='step', label="SM")

    print(val_of_bins_data1)

    ax1.legend()

    # Set ratio where val_of_bins_data2 is not zero
    ratio = np.divide(val_of_bins_data2,
                      val_of_bins_data1,
                      where=(val_of_bins_data1 != 0))

    print("ratio:", ratio)

    # Compute error on ratio (null if cannot be computed)
    error = np.divide(val_of_bins_data2 * np.sqrt(val_of_bins_data1) + val_of_bins_data1 * np.sqrt(val_of_bins_data2),
                      np.power(val_of_bins_data1, 2),
                      where=(val_of_bins_data1 != 0))

    print("error:", error)
 
    # Add the ratio on the existing plot
    # Add an histogram of errors
    ax3 = fig.add_subplot(3,1,3)
    ax3.set_ylabel('BSM/SM')
    ax3.set_xticks(np.arange(550,3550+300,300)) #set ticks

    bincenter = 0.5 * (edges_of_bins_data1[1:] + edges_of_bins_data1[:-1])
    ax3.errorbar(bincenter, ratio, yerr=error, fmt='o', color='k') # ratio with error
    # ax3.hist(error, nbins, range=(550,3550), histtype='step') #error 
    ax3.set_ylim(-0.5,3.)
 

#    ax3.set_xlabel('error')
#    ax3.set_ylabel('count')



#     main_ax_artists, subplot_ax_artists = hist_1.plot_ratio(
#         hist_2,
#         rp_ylabel="BSM/SM", # y axis label (lower plot) 
#         rp_num_label=thelabel,
#         rp_denom_label="SM",
#         rp_uncert_draw_type="line",  # line or bar
#         rp_uncertainty_type="poisson",
#         rp_ylim = [-0.5,5.5]
#     )
# #    subplot_ax_artists.set_ylabel("cross section [fb/100 GeV]") # y axis label (top plot)
    plt.savefig(f"{outdir}/{hist_name}.pdf")



def analyze(processo,oppe,valu):
    '''Event loop + histogram filling'''

#    lhe_file = os.path.join('/afs', 'cern.ch', 'user', 'a', 'acappati', 'work', 'ZZH', '220414_process1_nocuts', 'MG5_aMC_v2_7_3_py3', processo, 'Events', 'run_' + oppe + '_' + valu + '_cuts', 'unweighted_events.lhe') # process 1
    lhe_file = os.path.join('/afs', 'cern.ch', 'user', 'a', 'acappati', 'work', 'ZZH', '220414_process3_nocuts', 'MG5_aMC_v2_7_3_py3', processo, 'Events', 'run_' + oppe + '_' + valu + '_nocuts', 'unweighted_events.lhe') # process 3
    lhe_file_gz = lhe_file + '.gz'

    # check if gzipped file exists
    if os.path.isfile(lhe_file_gz):
        # open the gzip (read+byte mode) as file_in context
        # open the recipient unzipped file (write+byte mode) as file_out context
        with gzip.open(lhe_file_gz, 'rb') as file_in, open(lhe_file, 'wb') as file_out:
            # copy content of zipped file to recipient
            shutil.copyfileobj(file_in, file_out)

    # check if unzipped file exists
    print('opening file ', lhe_file)
    if os.path.isfile(lhe_file):
        reader = LHEReader(lhe_file)
    else:
        # throw error
        raise FileNotFoundError(f'{lhe_file} not found!')
    
    # array delle masse
    mass_array = []
 
    # loop over events 
    for event in reader:
        # Find bosons
        tops = filter(
            lambda x: abs(x.pdgid) in (23,24,25),
            event.particles
        )

        # Sum over top four-momenta in the event
        combined_p4 = None
        for p4 in map(lambda x: x.p4(), tops):
            if combined_p4:
                combined_p4 += p4
            else:
                combined_p4 = p4

        mass_array.append(combined_p4.mass) # fill mass array

    return mass_array


def main():

    if (len(sys.argv) < 3):
        print("specify the process, operator and valu please")
        sys.exit(1)

    # --- out directory
    processo = sys.argv[1]
    oppe = sys.argv[2]
    valu = sys.argv[3]
    valu2 = sys.argv[4]
    xs = sys.argv[5]
    xs2 = sys.argv[6]
    label = sys.argv[7]

    #histograms = analyze('/afs/cern.ch/work/c/covarell/mg5_amcatnlo/test-dim8-zzh/MG5_aMC_v2_7_3_py3/vbf-hh-mhhcut/Events/run_05/unweighted_events.lhe')
    #histograms = analyze('/afs/cern.ch/user/c/covarell/work/mg5_amcatnlo/dim8-hh/MG5_aMC_v2_7_3_py3/vbf-wpmz-4f/Events/run_FM4_20_cutshistat/unweighted_events.lhe')
    data1 = analyze(processo,oppe,valu)
    data2 = analyze(processo,oppe,valu2)
    plot(data1,data2,oppe,valu2,xs,xs2,label)

if __name__=="__main__":
    main()
    exit(0)
