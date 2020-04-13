
date = '180221'; import os
todays_folder = 'Y:\Lab\MOT\Daily\\'+date[:4]+'\\'+date;
os.chdir(todays_folder+'\\scripts\\'+'analysis ipython')
get_ipython().magic(u'run standard_functions.ipynb')
from scipy.stats import poisson


def plotbayesian(
    path, conditions,
    filter_recap_q      = False,
    filter_no_recap_q   = False,
    limit_file_q        = False,
    limitNumFiles       = 1,
    filter_ai_q         = False,
    lock_cutoff         = 700,
    lock_cutoff_lower   = 0,
    filter_coupling_q   = False,
    time_beg_coupling   = 0,
    time_end_coupling   = 100000,
    count_threshold     = 10,
    N_bins              = 250,
    hist_range_end      = 500*2500,
    disable_lines_q     = False,
    means               = None
    ):
    #----------- default values that don't change often and initialize file
    os.chdir(path)
    h5filename      = path+'.h5'
    f               = h5.File(h5filename,'r')
    # ----------get files to be used
    timesFiles      = [filename for filename in os.listdir(path) if filename.startswith('times_')]
    paramFiles      = [filename for filename in os.listdir(path) if filename.startswith('params_')]
    aiFiles      = [filename for filename in os.listdir(path) if filename.startswith('ai_')]
    #--------------------get lock data
    ai_data = {}
    for aifile in aiFiles:
        midpoint, ai_trace, data = aitrace(aifile,AI_startPt,AI_endPt)
        ai_data[aifile.split(",")[1][4:]]=(midpoint, ai_trace, data)
    # -----------seq filter if there is more than one variable

    print "Conditions for the selected point: ",conditions
    selected = filterSeq(conditions,path,paramFiles,limit_file_q,limitNumFiles)
    print "Number of sequences corresponding to the condition: ",len(selected)

    # ------------- initialize counters
    countNumSeqDone = 0;
    total_N_atoms = 0
    data_atoms = []
    #-------- loop through sequences
    for seq in selected:

        current_seq_name = seq[3:]

        paramFile = getParamFile(seq,paramFiles)                  # get the param file
        # if there is no param file we can't do anything, so skip straight away to next seq
        if paramFile==False:             # no match
            print 'Skipping  seq ', current_seq_name," :no param file"
            continue
        try:
            time, aid, pid = getIDsFromSeq(seq,channel,f)         # get time, aid, pid
        except KeyError:
            print 'Skipping  seq ', current_seq_name," :h5 error"
            continue

        #-----------------convert time to ns using picoharp's resolution
        time = np.multiply(time,0.004)
        #-------------now we have aid, pid, times, apply filters
        if (filter_recap_q):

            timesFile = getTimesFile(seq,timesFiles)
            if timesFile == False: # no match
                print 'Skipping  seq ', current_seq_name," :no times file"
                continue
            time, aid, pid = filter_recapture(path,timesFile,paramFile,time,aid,pid)
        if (filter_no_recap_q):

            timesFile = getTimesFile(seq,timesFiles)
            if timesFile == False: # no match
                print 'Skipping  seq ', current_seq_name," :no times file"
                continue
            time, aid, pid = filter_no_recapture(path,timesFile,paramFile,time,aid,pid)

        if (filter_ai_q):
            try:
                current_lock = ai_data[current_seq_name][1]
                time = [time[i] for i in np.arange(len(time)) if (current_lock[aid[i]-1]<lock_cutoff)and(lock_cutoff_lower<current_lock[aid[i]-1])]
                aid_new  = [aid[i] for i in np.arange(len(aid)) if (current_lock[aid[i]-1]<lock_cutoff)and(lock_cutoff_lower<current_lock[aid[i]-1])]
                pid  = [pid[i] for i in np.arange(len(pid)) if (current_lock[aid[i]-1]<lock_cutoff)and(lock_cutoff_lower<current_lock[aid[i]-1])]
                aid = aid_new
            except KeyError:
                print 'Skipping  seq ', current_seq_name," :no ai file"
                continue
        if (filter_coupling_q):
            time, aid, pid = filter_coupling(time_beg_coupling,time_end_coupling,count_threshold,time,aid,pid)

# ----------add total atom number for normalization to this sequence
        total_N_atoms=total_N_atoms+len(set(aid))
# ---------separate data per into atoms ---------------------------
        for atom in set(aid):
        # which pulses correspond to this atom id
            indexes = [i for i, j in enumerate(aid) if j == atom]
            #print 'indexes ',indexes
            pulseTimes = [time[i] for i in indexes]
            # if we want to filter out atoms that had more than filter_coupling_threshold counts overall
            data_atoms.append(np.histogram(pulseTimes,bins=N_bins,range = (0,hist_range_end))[0])

        # -----------count how many sequences were done
        countNumSeqDone = countNumSeqDone+1
        if limit_file_q:
            if countNumSeqDone>limitNumFiles:
                break
    f.close()
    #---------------------------------------- start with plotting ---------------------------
    print "Number of sequences: ",countNumSeqDone
    print "Total number of atoms: ",total_N_atoms

    mean_counts = np.mean(flatten(data_atoms))

    #colors= itertools.cycle(["#FF4848","#003F87","black","yellow","green","magenta","chocolate","orange"])
    colors= itertools.cycle(["#FF4848","#003F87","black"])
    x_axis = np.linspace(0,hist_range_end,N_bins)
    x_axis = np.divide(x_axis,1000.0) # in mus
    bin_size = hist_range_end/(N_bins*1000.0) # in mus


    bins_for_fit=25
    c = flatten(data_atoms)
    points,bin_edges = np.histogram(c,bins_for_fit)

    x_axis2 = [(bin_edges[i]+bin_edges[i+1])/2.0 for i in range(len(bin_edges)-1)]

    if means is None:
        p0 = [np.mean(c),np.std(c),np.mean(c),np.std(c),1,1]
    else:
        p0 = [means[0],np.sqrt(means[0]),means[1],np.sqrt(means[1]),1,1]
    try:

        popt,pcov = curve_fit(double_gauss,x_axis2,points,p0=p0)

        print "Optimal fit parameters "
        print "[mean,std,shot_noise_std/std,a] = ", popt[0],np.abs(popt[1]),np.sqrt(popt[0])/np.abs(popt[1]),popt[4]
        print "[mean,std,shot_noise_std/std,a] = ", popt[2],np.abs(popt[3]),np.sqrt(popt[2])/np.abs(popt[3]),popt[5]

    except RuntimeError:
        print "Can't find optimal parameters for double Gaussian fit"


    print "Plotting individual traces overlapped with losing probability"
    for jj in range(len(data_atoms)):#[range(len(data_atoms))[-1]]: #13
        data_by_atoms = data_atoms[jj]
        r1 =   popt[0]/float(bin_size)
        r2 = popt[2]/float(bin_size)
        # r1 should be the higher rate
        if r1<r2:
            helper=r1
            r1=r2
            r2=helper

        Bayesian_losing(x_axis,data_by_atoms,r1,r2)

def adjust_spines(ax,spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward',10)) # outward by 10 points
        else:
            spine.set_color('none') # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])



import matplotlib.ticker as ticker
def Bayesian_losing(time,counts,r1,r2):
    bin_size = time[-1]/len(time)
    T = time[-1]
    prob_lost=[]
    for i in range(len(time)):
        t_lost=bin_size*(i+1)
        N_photons_before = np.sum(counts[0:i+1])
        N_photons_after = np.sum(counts[i:])
        expected_N_photons_before = r1*t_lost
        expected_N_photons_after = r2*(T-t_lost)

        p1 = normal_prob(N_photons_before,expected_N_photons_before,np.sqrt(expected_N_photons_before))
        p2 = normal_prob(N_photons_after,expected_N_photons_after,np.sqrt(expected_N_photons_after))
        prob = p1*p2
        prob_lost.append(prob)

    fig = plt.figure(figsize=(10,6))

    ax = fig.add_subplot(2,1,1)
    ax.plot(time,counts,color="r",marker='o',linestyle='--')
    ax.set_ylabel("counts per bin")
    ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    adjust_spines(ax,['left'])

    ax = fig.add_subplot(2,1,2)
    normfactor = 1/np.sum(prob_lost[:-1])
    prob_lost = np.multiply(prob_lost,normfactor)

    ax.plot(time,prob_lost,color = "#003F87",marker='o',linestyle='-.')
    ax.yaxis.set_major_locator(plt.MaxNLocator(4))

    ax.set_ylabel("prob. of losing atom")
    ax.set_xlabel("time [us]")
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    adjust_spines(ax,['left','bottom'])

    plt.savefig(todays_folder+"\\bayesian3")
    plt.show()


path = r'Y:\Lab\MOT\Daily\1710\171020\pc12';
xName = "blaster_detuning_MHz"; plotBy = "dump_atom_q";
lockcut = 450; t_end = 50000; couplingcut = 60;

for z in ["70.0"]:
    conditions1 = [(xName,z),(plotBy,"0")]

    plotbayesian(
        path, conditions1,
        filter_recap_q    = False,
        filter_no_recap_q = False,
        limit_file_q    = False,
        limitNumFiles   = 1,
        filter_ai_q=True,
        lock_cutoff = lockcut,lock_cutoff_lower = 0,
        filter_coupling_q = True,
        time_beg_coupling=0,time_end_coupling=t_end,count_threshold=couplingcut,
        N_bins = 40,hist_range_end= 800000,disable_lines_q=True,means=[20,65]
    )
