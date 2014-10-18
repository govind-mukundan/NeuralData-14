#
#  NAME
#    problem_set1.py
#
#  DESCRIPTION
#    Open, view, and analyze raw extracellular data
#    In Problem Set 1, you will write create and test your own spike detector.
#

import numpy as np
import matplotlib.pylab as plt

def load_data(filename):
    """
    load_data takes the file name and reads in the data.  It returns two 
    arrays of data, the first containing the time stamps for when they data
    were recorded (in units of seconds), and the second containing the 
    corresponding voltages recorded (in units of microvolts - uV)
    """
    data = np.load(filename)[()];
    return np.array(data['time']), np.array(data['voltage'])
    
def bad_AP_finder(time,voltage):
    """
    This function takes the following input:
        time - vector where each element is a time in seconds
        voltage - vector where each element is a voltage at a different time
        
        We are assuming that the two vectors are in correspondance (meaning
        that at a given index, the time in one corresponds to the voltage in
        the other). The vectors must be the same size or the code
        won't run
    
    This function returns the following output:
        APTimes - all the times where a spike (action potential) was detected
         
    This function is bad at detecting spikes!!! 
        But it's formated to get you started!
    """
    
    #Let's make sure the input looks at least reasonable
    if (len(voltage) != len(time)):
        print "Can't run - the vectors aren't the same length!"
        APTimes = []
        return APTimes
    
    numAPs = np.random.randint(0,len(time))//10000 #and this is why it's bad!!
    print "Number of APs found = " + str(numAPs)
 
    # Now just pick 'numAPs' random indices between 0 and len(time)
    APindices = np.random.randint(0,len(time),numAPs)
    
    # By indexing the time array with these indices, we select those times
    APTimes = time[APindices]
    
    # Sort the times
    APTimes = np.sort(APTimes)
    
    return APTimes
    
def good_AP_finder(time,voltage):
    """
    This function takes the following input:
        time - vector where each element is a time in seconds
        voltage - vector where each element is a voltage at a different time
        
        We are assuming that the two vectors are in correspondance (meaning
        that at a given index, the time in one corresponds to the voltage in
        the other). The vectors must be the same size or the code
        won't run
    
    This function returns the following output:
        APTimes - all the times where a spike (action potential) was detected
    """
    APTimes = []
    #Let's make sure the input looks at least reasonable
    if (len(voltage) != len(time)):
        print "Can't run - the vectors aren't the same length!"
        return APTimes
    
    ##Your Code Here!
    Window_ms = 1
    Sample_Rate_Hz = 1/(time[1] - time[0])
    Window_cnt = int(Sample_Rate_Hz * Window_ms/1000 + 0.5)
    voltage_rect = np.abs(voltage)
    Threshold = max(voltage_rect)/2 #half the max voltage is the threshold
    print "Sampling Rate = " +str(Sample_Rate_Hz)
    print "Threshold = " + str(Threshold)
    print "Samples in the Window = " + str(Window_cnt)
    voltage_mw = np.ones(len(voltage))
    
    # Find the moving average of the signal over the window
#    for j in range(Window_cnt,len(voltage)): # Initial samples are 1
#        for i in range(0,Window_cnt):
#            voltage_mw[j] = voltage_mw[j] + voltage_rect[j - i]
#    
#    voltage_mw = voltage_mw/Window_cnt
    #plt.plot(time,voltage_mw)
    
    # Find peaks above the threshold in the averaged signal
    # A peak is detected if the signal goes above and below the threshold within the window interval
    LookForFall = False
    Peak_Indices = []
    Pk_Interval = Window_cnt
    for j in range(0,len(voltage_rect)):
        Pk_Interval = Pk_Interval + 1
        if (LookForFall == False):
            if(Pk_Interval > Window_cnt and voltage_rect[j] > Threshold):
                LookForFall = True
                Pk_Interval = 0
        else:
            if(Pk_Interval > Window_cnt): # Timeout
                LookForFall = False
                print "Timeout, couldnt find falling edge"
            elif(voltage_rect[j] < Threshold): #Falling edge found
                LookForFall = False
                APTimes.append(time[0] + j/Sample_Rate_Hz) #Convert sample number into time 
                print "Peaks Found at " + str(j) + " = " + str(time[0] + j/Sample_Rate_Hz)
    
    print "Total Peaks = " + str(len(APTimes))
        
       
    return APTimes
    

def get_actual_times(dataset):
    """
    Load answers from dataset
    This function takes the following input:
        dataset - name of the dataset to get answers for

    This function returns the following output:
        APTimes - spike times
    """    
    return np.load(dataset)
    
def detector_tester(APTimes, actualTimes):
    """
    returns percentTrueSpikes (% correct detected) and falseSpikeRate
    (extra APs per second of data)
    compares actual spikes times with detected spike times
    This only works if we give you the answers!
    """
    
    JITTER = 0.025 #2 ms of jitter allowed
    
    #first match the two sets of spike times. Anything within JITTER_MS
    #is considered a match (but only one per time frame!)
    
    #order the lists
    detected = np.sort(APTimes)
    actual = np.sort(actualTimes)
    
    #remove spikes with the same times (these are false APs)
    temp = np.append(detected, -1)
    detected = detected[plt.find(plt.diff(temp) != 0)]
 
    #find matching action potentials and mark as matched (trueDetects)
    trueDetects = [];
    for sp in actual:
        z = plt.find((detected >= sp-JITTER) & (detected <= sp+JITTER))
        if len(z)>0:
            for i in z:
                zz = plt.find(trueDetects == detected[i])
                if len(zz) == 0:
                    trueDetects = np.append(trueDetects, detected[i])
                    break;
    percentTrueSpikes = 100.0*len(trueDetects)/len(actualTimes)
    
    #everything else is a false alarm
    totalTime = (actual[len(actual)-1]-actual[0])
    falseSpikeRate = (len(APTimes) - len(actualTimes))/totalTime
    
    print 'Action Potential Detector Performance performance: '
    print '     Correct number of action potentials = ' + str(len(actualTimes))
    print '     Percent True Spikes = ' + str(percentTrueSpikes)
    print '     False Spike Rate = ' + str(falseSpikeRate) + ' spikes/s'
    print 
    return {'Percent True Spikes':percentTrueSpikes, 'False Spike Rate':falseSpikeRate}
    
    
def plot_spikes(time,voltage,APTimes,titlestr):
    """
    plot_spikes takes four arguments - the recording time array, the voltage
    array, the time of the detected action potentials, and the title of your
    plot.  The function creates a labeled plot showing the raw voltage signal
    and indicating the location of detected spikes with red tick marks (|)
    """
    plt.figure()
    
    ##Your Code Here  
    plt.plot(time,voltage)
    plt.hold = True
    plt.xlabel("Time (S)")
    plt.ylabel("Voltage (uV)")
    plt.title(titlestr)
    # Generate an impulse array with impulses at APPoints
    APv = np.ones(len(APTimes))
    APv = APv * (max(voltage) + 10)
    plt.plot(APTimes, APv, 'r|')  # plot x and y using blue circle markers
    
    plt.show()
    
def plot_waveforms(time,voltage,APTimes,titlestr):
    """
    plot_waveforms takes four arguments - the recording time array, the voltage
    array, the time of the detected action potentials, and the title of your
    plot.  The function creates a labeled plot showing the waveforms for each
    detected action potential
    """
   
    plt.figure()

    ## Your Code Here   
    plt.hold = True
    plt.ylabel("Voltage (uV)")
    plt.xlabel("Time (S)")
    plt.title(titlestr)
    
    APT = []
    Sample_Rate_Hz = 1/(time[1] - time[0])
    # Plot the APs detected over a time period of +/- 3mS
    Sample_Cnt = int(Sample_Rate_Hz * 6 / 1000)
    Time = range(0,Sample_Cnt)
    # Convert sample numbers into time starting with -3mS
    Time[:] = [-0.3 + x/Sample_Rate_Hz for x in Time]
    print "sample count = " + str(len(Time))
    #APTimes = APTimes * Sample_Rate_Hz #convert into index
    APT[:] = [((x - time[0])*Sample_Rate_Hz) for x in APTimes] 
    print "Len AP = " + str(len(APTimes))
    for i in range(0,len(APTimes)):
       # Get the data to plot
       data = voltage[APT[i] - Sample_Cnt/2 : APT[i] + Sample_Cnt/2]
       #print(len(data))
       if(len(Time) == len(data)):
           plt.plot(Time,data,'b')
       #
   
    plt.show()
    

        
##########################
#You can put the code that calls the above functions down here    
if __name__ == "__main__":
#    t,v = load_data('spikes_example.npy')    
#    actualTimes = get_actual_times('spikes_example_answers.npy')
    t,v = load_data('spikes_challenge.npy')    
    actualTimes = get_actual_times('spikes_hard_practice_answers.npy')
    APTime = good_AP_finder(t,v)
    plot_spikes(t,v,APTime,'Action Potentials in Raw Signal')
    plot_waveforms(t,v,APTime,'Waveforms')
    detector_tester(APTime,actualTimes)


