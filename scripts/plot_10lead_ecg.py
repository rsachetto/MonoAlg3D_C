# ===============================================================================
# Authors: Jenny Wang, Julia Camps and Lucas Berg
# Version: 24/06/24
# ===============================================================================
# This is a stand alone script that plots ECGs from body surface potentials
# simulated from monodomain simulations.
# The normalisation strategy is the same as in ecg_functions.py.
# ===============================================================================
# Basic script to plot the full ECG.
# ===============================================================================
import sys
import numpy as np
import scipy
from scipy import signal
#from scipy import stats
import matplotlib
from matplotlib import pyplot as plt
# Helper functions from ecg_functions.py:

# Global variables
global filtering
filtering = True
global normalise
normalise = True
global zero_align
zero_align = True
global frequency
frequency = 1000
global high_freq_cut
high_freq_cut = 150
global low_freq_cut
low_freq_cut = 0.5
global max_len_ecg
max_len_ecg = 450
global max_len_qrs
max_len_qrs = 200
global nb_leads
nb_leads = 8
global qrs_onset
qrs_onset = 0
#global casename
#casename = 'DTI004'
global monoalg_activation_offset
monoalg_activation_offset=29.43 # DTI004
global ref_ecg_offset
ref_ecg_offset = 101 # DTI004

def filter_butterworth_ecg(b, a, ecg):
    return signal.filtfilt(b, a, ecg)

def __filter_ecg(original_ecg, low_a_filtfilt, low_b_filtfilt, high_a_filtfilt, high_b_filtfilt):
    # First we filter out the low frequencies using a high-pass filter and the lower thresholds
    processed_ecg = filter_butterworth_ecg(b=low_b_filtfilt, a=low_a_filtfilt, ecg=original_ecg)
    # Secondly we filter out the high frequencies using a low-pass filter and the higher thresholds
    return filter_butterworth_ecg(b=high_b_filtfilt, a=high_a_filtfilt, ecg=processed_ecg)

def zero_align_ecg(original_ecg):
    return original_ecg - (original_ecg[:, 0:1] + original_ecg[:, -2:-1]) / 2  # align at zero


def __preprocess_ecg_without_normalise(original_ecg, filtering, zero_align, frequency, high_freq_cut, low_freq_cut):
    processed_ecg = original_ecg
    if filtering:
        high_w = high_freq_cut / (frequency / 2)  # Normalize the frequency
        high_b_filtfilt, high_a_filtfilt = signal.butter(4, high_w,'low')  # Butterworth filter of fourth order.
        low_w = low_freq_cut / (frequency / 2)  # Normalize the frequency
        low_b_filtfilt, low_a_filtfilt = signal.butter(4, low_w, 'high')
        processed_ecg = __filter_ecg(processed_ecg, low_a_filtfilt=low_a_filtfilt, low_b_filtfilt=low_b_filtfilt,
                                     high_a_filtfilt=high_a_filtfilt, high_b_filtfilt=high_b_filtfilt)
    if zero_align:
        processed_ecg = zero_align_ecg(processed_ecg)
    return processed_ecg

def __set_reference_lead_is_positive(max_qrs_end, reference_ecg):
    approx_qrs_end = min(reference_ecg.shape[1], max_qrs_end)  # Approximate end of QRS.
    reference_lead_max = np.absolute(np.amax(reference_ecg[:, :approx_qrs_end], axis=1))
    reference_lead_min = np.absolute(np.amin(reference_ecg[:, :approx_qrs_end], axis=1))
    reference_lead_is_positive = reference_lead_max >= reference_lead_min
    reference_amplitudes = np.zeros(shape=nb_leads, dtype=np.float64)
    reference_amplitudes[reference_lead_is_positive] = reference_lead_max[reference_lead_is_positive]
    reference_amplitudes[np.logical_not(reference_lead_is_positive)] = reference_lead_min[
        np.logical_not(reference_lead_is_positive)]
    # if verbose:
    #     print('reference_lead_is_positive')
    #     print(reference_lead_is_positive)
    #     print('reference_amplitudes')
    #     print(reference_amplitudes)
    return reference_lead_is_positive  # Have some R progression by normalising by the
    # largest absolute amplitude lead

def __normalise_ecg_based_on_rwave_8_leads(original_ecg, qrs_onset, max_len_qrs):
    if nb_leads != 8 or original_ecg.shape[0] != 8:
        raise(Exception, 'This function is hardcoded for the specific ECG configuration: I, II, V1, ..., V6')
    # print('Normalising ECG ', original_ecg.shape)
    #approx_qrs_end = min(reference_ecg.shape[1], max_len_qrs+qrs_onset)  # Approximate end of QRS.
    # approx_qrs_width = min(original_ecg.shape[1], max_len_qrs)  # This approximation is more robust.
    # print('approx_qrs_end ', approx_qrs_end)
    # print(np.amax(original_ecg[:, :approx_qrs_end]))
    #reference_amplitudes = np.empty(shape=nb_leads, dtype=np.float64)
    #reference_amplitudes[reference_lead_is_positive] = np.absolute(
    #    np.amax(original_ecg[:, :approx_qrs_end], axis=1)[
    #        reference_lead_is_positive])
    #reference_amplitudes[np.logical_not(reference_lead_is_positive)] = \
    #    np.absolute(np.amin(original_ecg[:, :approx_qrs_end], axis=1))[np.logical_not(
    #        reference_lead_is_positive)]
    normalised_ecg = np.zeros(original_ecg.shape)
    normalised_ecg[:2, :] = original_ecg[:2, :]
    normalised_ecg[2:nb_leads, :] = original_ecg[2:nb_leads, :]
    return normalised_ecg

def preprocess_ecg(original_ecg, filtering, zero_align, frequency, high_freq_cut, low_freq_cut, max_len_qrs, qrs_onset):
    processed_ecg = original_ecg
    if filtering:
        high_w = high_freq_cut / (frequency / 2)  # Normalize the frequency
        high_b_filtfilt, high_a_filtfilt = signal.butter(4, high_w,'low')  # Butterworth filter of fourth order.
        low_w = low_freq_cut / (frequency / 2)  # Normalize the frequency
        low_b_filtfilt, low_a_filtfilt = signal.butter(4, low_w, 'high')
        processed_ecg = __filter_ecg(processed_ecg, low_a_filtfilt=low_a_filtfilt, low_b_filtfilt=low_b_filtfilt,
                                     high_a_filtfilt=high_a_filtfilt, high_b_filtfilt=high_b_filtfilt)
    if zero_align:
        processed_ecg = zero_align_ecg(processed_ecg)
    if normalise:
        processed_ecg = __normalise_ecg_based_on_rwave_8_leads(original_ecg=processed_ecg, qrs_onset=qrs_onset, max_len_qrs=max_len_qrs)
    return processed_ecg

def import_simulated_ecg_8leads_raw(filename, monoalg_activation_offset):
    data = np.genfromtxt(filename, skip_footer=1)
    #t = data[:, 0] - monoalg_activation_offset
    t = data[:, 0]
    LA = data[:, 1]
    RA = data[:, 2]
    LL = data[:, 3]
    RL = data[:, 4]
    V1 = data[:, 5]
    V2 = data[:, 6]
    V3 = data[:, 7]
    V4 = data[:, 8]
    V5 = data[:, 9]
    V6 = data[:, 10]

    # Ealuate Wilson's central terminal
    VW = 1.0 / 3.0 * (RA + LA + LL)

    # Evaluate simulated ECG lead traces
    V1 = V1 - VW
    V2 = V2 - VW
    V3 = V3 - VW
    V4 = V4 - VW
    V5 = V5 - VW
    V6 = V6 - VW
    I = LA - RA
    II = LL - RA
    III = LL - LA
    aVL = LA - (RA + LL) / 2.0
    aVF = LL - (LA + RA) / 2.0
    aVR = RA - (LA + LL) / 2.0
    ecgs = np.vstack((I, II, V1, V2, V3, V4, V5, V6))
    return t, ecgs

def visualise_ecgs(nb_leads, simulated_ecgs, simulated_t,  lead_names, casename, output_filename):
    nb_cols = (nb_leads * 2) ** 0.5
    if nb_cols - int(nb_cols) == 0. and nb_cols / 2 - int(nb_cols / 2) == 0.:
        nb_rows = nb_cols / 2
    else:
        # Try to make 2 rows and the necessary columns
        nb_cols = nb_leads / 2
        if nb_cols - int(nb_cols) == 0. and nb_cols / 2 - int(nb_cols / 2) == 0.:
            nb_rows = nb_cols / 2
        else:
            nb_cols = nb_leads
            nb_rows = 1
    fig, axes = plt.subplots(int(nb_rows), int(nb_cols), figsize=(20, 10))
    axes = np.reshape(axes, nb_leads)
    for lead_i in range(nb_leads):
        time_steps = np.arange(simulated_ecgs.shape[1])
        #axes[lead_i].plot(time_steps[reference_ecg_offset:]-reference_ecg_offset, reference_ecg[lead_i, reference_ecg_offset:], label='Clinical', color='orangered', linewidth=3., linestyle='--')
        #axes[lead_i].plot(simulated_t, simulated_ecgs[lead_i, :], color='k', label='Simulation', linewidth=2.)
        axes[lead_i].plot(simulated_t, simulated_ecgs[lead_i,:], color='k', label='Simulated', linewidth=2.)
        axes[lead_i].set_title(lead_names[lead_i], fontsize=20)
        axes[lead_i].set_ylim([-0.1, 0.1])
        for tick in axes[lead_i].xaxis.get_major_ticks():
            tick.label1.set_fontsize(12)
        for tick in axes[lead_i].yaxis.get_major_ticks():
            tick.label1.set_fontsize(12)
        for axis in ['top','bottom','left','right']:
            axes[lead_i].spines[axis].set_linewidth(2)
    axes[nb_leads-1].legend(loc='upper right', fontsize=12)
    fig.suptitle("Simulated ECG traces", fontsize=24)
    plt.show()
    #plt.savefig(output_filename, dpi=300)

def calculate_pearson_correlation(nb_leads, reference_ecg, simulated_ecgs, simulated_t,  lead_names, casename, output_filename, reference_ecg_offset):
    
    ref_time_steps = np.arange(reference_ecg.shape[1])
    ref_time_steps = ref_time_steps[reference_ecg_offset:672]-reference_ecg_offset

    pearson_corr_arr = []
    for lead_i in range(nb_leads):
        ref_ecg = reference_ecg[lead_i,reference_ecg_offset:672]
        #print("%d %d" % (len(ref_time_steps), len(ref_ecg)))

        aprox_time_steps = []
        for i in range(len(simulated_t)):
            if (i % 10 == 0):
                aprox_time_steps.append(simulated_t[i])
        aprox_time_steps = np.array(aprox_time_steps)
        aprox_ecg = []
        for i in range(len(simulated_ecgs[lead_i])):
            if (i % 10 == 0):
                aprox_ecg.append(simulated_ecgs[lead_i,i])
        
        aprox_time_steps = aprox_time_steps[int(monoalg_activation_offset):]-int(monoalg_activation_offset)
        aprox_ecg = aprox_ecg[int(monoalg_activation_offset):]

        pearson_corr = scipy.stats.pearsonr(ref_ecg,aprox_ecg)
        pearson_corr_arr.append(pearson_corr[0])
    
    return np.array(pearson_corr_arr)

def main ():
    if (len(sys.argv) != 4):
        print("="*100)
        print("Usage:> python %s <meshname> <input_simulation_ecg> <output_png_ecg>" % (sys.argv[0]))
        print("="*100)
        sys.exit(1)

    ###################################################################################################################
    ## Main script
    casename = sys.argv[1]
    monoalg_ecg_filename = sys.argv[2]
    output_ecg_filename = sys.argv[3]

    # Preprocess the clinical data
    #print('Preprocessing clinical ECGs')
    #clinical_data_filename_path = "inputs/%s_clinical_full_ecg.csv" % (casename)
    #reference_ecg = np.genfromtxt(clinical_data_filename_path, delimiter=',')
    #reference_ecg = __preprocess_ecg_without_normalise(original_ecg=reference_ecg, filtering=filtering, zero_align=zero_align,
    #                                                frequency=frequency, low_freq_cut=low_freq_cut, high_freq_cut=high_freq_cut)
    #max_qrs_end_aux = qrs_onset + max_len_qrs
    #reference_lead_is_positive = __set_reference_lead_is_positive(max_qrs_end=max_qrs_end_aux, reference_ecg=reference_ecg)
    #reference_ecg = preprocess_ecg(original_ecg=reference_ecg, reference_ecg=reference_ecg, filtering=filtering, zero_align=zero_align,
    #                                frequency=frequency, low_freq_cut=low_freq_cut, high_freq_cut=high_freq_cut, max_len_qrs=max_len_qrs,
    #                            qrs_onset=qrs_onset, reference_lead_is_positive=reference_lead_is_positive)

    # Read in simulated ECGs
    #print('Importing and preprocessing simulated monoAlg3D ECGs')
    simulated_t, simulated_ecgs_8leads = import_simulated_ecg_8leads_raw(filename=monoalg_ecg_filename, monoalg_activation_offset=monoalg_activation_offset)
    #simulated_t, simulated_ecgs_8leads = import_simulated_ecg_8leads_raw(filename=casename+'_ecg.txt', monoalg_activation_offset=37)
    simulation_frequency = simulated_ecgs_8leads.shape[1]
    max_len_qrs_aux = 200 * 4
    processed_simulated_ecgs = preprocess_ecg(original_ecg=simulated_ecgs_8leads, 
                                            filtering=filtering, zero_align=zero_align, frequency=simulation_frequency, low_freq_cut=low_freq_cut, high_freq_cut=high_freq_cut, max_len_qrs=max_len_qrs_aux,
                                            qrs_onset=qrs_onset)

    # Plot together
    print('Visualising ECG ...')
    visualise_ecgs(nb_leads=nb_leads, simulated_ecgs=processed_simulated_ecgs, simulated_t=simulated_t,
                lead_names=['I', 'II', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'], casename=casename, output_filename=output_ecg_filename)

if (__name__ == "__main__"):
    main()
