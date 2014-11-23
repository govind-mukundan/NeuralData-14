#
#  NAME
#    problem_set4.py
#
#  DESCRIPTION
#    In Problem Set 4, you will classify EEG data into NREM sleep stages and
#    create spectrograms and hypnograms.
#
from __future__ import division
import numpy as np
import matplotlib.pylab as plt
import matplotlib.mlab as m


def load_examples(filename):
    """
    load_examples takes the file name and reads in the data.  It returns an
    array containing the 4 examples of the 4 stages in its rows (row 0 = REM;
    1 = stage 1 NREM; 2 = stage 2; 3 = stage 3 and 4) and the sampling rate for
    the data in Hz (samples per second).
    """
    data = np.load(filename)
    return data['examples'], int(data['srate'])

def load_eeg(filename):
    """
    load_eeg takes the file name and reads in the data.  It returns an
    array containing EEG data and the sampling rate for
    the data in Hz (samples per second).
    """
    data = np.load(filename)
    return data['eeg'], int(data['srate'])

def load_stages(filename):
    """
    load_stages takes the file name and reads in the stages data.  It returns an
    array containing the correct stages (one for each 30s epoch)
    """
    data = np.load(filename)
    return data['stages']

def plot_example_psds(example,rate):
    """
    This function creates a figure with 4 lines to show the overall psd for 
    the four sleep examples. (Recall row 0 is REM, rows 1-3 are NREM stages 1,
    2 and 3/4)
    """
    plt.figure()
    
    ##YOUR CODE HERE    
    
    return

def plot_example_spectrograms(example,rate):
    """
    This function creates a figure with spectrogram sublpots to of the four
    sleep examples. (Recall row 0 is REM, rows 1-3 are NREM stages 1,
    2 and 3/4)
    """
    plt.figure()
    
    ###YOUR CODE HERE
    y_lim = 40
    plt.title('Spectrogram')
    bin_space = 512 #30*rate # A typical window size is 30 seconds
    plt.subplot(411)
    plt.specgram(examples[0]/np.sum(examples[0]),NFFT=bin_space,Fs=srate)
    plt.ylim((0,y_lim))
    plt.title ('REM')
    plt.subplot(412)
    plt.title ('Stage 1 NREM')
    plt.specgram(examples[1]/np.sum(examples[1]),NFFT=bin_space,Fs=srate)
    plt.ylim((0,y_lim))
    plt.subplot(413)
    plt.title ('Stage 2 NREM')
    plt.specgram(examples[2]/np.sum(examples[2]),NFFT=bin_space,Fs=srate)
    plt.ylim((0,y_lim))
    plt.subplot(414)
    plt.title ('Stage 3/4 NREM')
    plt.specgram(examples[3]/np.sum(examples[3]),NFFT=bin_space,Fs=srate)
    plt.ylim((0,y_lim))
    plt.show();
    
    return
      
            
def classify_epoch(epoch,rate):
    """
    This function returns a sleep stage classification (integers: 1 for NREM
    stage 1, 2 for NREM stage 2, and 3 for NREM stage 3/4) given an epoch of 
    EEG and a sampling rate.
    """

    ###YOUR CODE HERE
    # Find the frequency content in the ranges - [2-7]Hz, and 10-15Hz and use the range [5-10]Hz as reference
    # If the 
    # 1. find the mag spectrum - psd
    print "---------"
    Pxx,freq = plt.psd(epoch/np.sum(epoch),NFFT=512,Fs=rate);
    Pxx = Pxx/np.sum(Pxx)
    psd_lf = np.max(10*np.log10(Pxx[2:7])) - np.min(10*np.log10(Pxx[2:7]));
    psd_ref = np.max(10*np.log10(Pxx[7:10])) - np.min(10*np.log10(Pxx[7:10]))
    psd_hf = np.max(10*np.log10(Pxx[12:14])) - np.min(10*np.log10(Pxx[12:14]))
    # Avg psd in low feq rang
    print "Avg 2-7Hz = " + str(psd_lf)
    print "Avg 7-10Hz = " + str(psd_ref)
    print "Avg 10-13Hz = " + str(psd_hf)
    
    alpha_pow_ratio = np.sum(Pxx[8:12]) /np.sum(Pxx[1:40]) 
    lf_pow_ratio = np.sum(Pxx[2:5]) /np.sum(Pxx[1:40]) 
    hf_pow_ratio = np.sum(Pxx[12:15]) /np.sum(Pxx[1:40])
    ref_pow_ratio = np.sum(Pxx[18:25]) /np.sum(Pxx[1:40])
    print "Low freq power = " + str(lf_pow_ratio * 100)
    print "High freq power = " + str(hf_pow_ratio * 100)
    print "Ref freq power = " + str(ref_pow_ratio * 100)
    print "Alpha freq power = " + str(alpha_pow_ratio * 100)
    
#    plt.figure();
#    plt.hist(Pxx[1:40]) # ,cumulative=True,histtype='step'
    
    pkh_hf1 = find_peak(Pxx[7:15])
    pkh_hf2 = find_peak(Pxx[0:5])
    pkh_ref = find_peak(Pxx[4:10])
    stage = 0
    if pkh_hf1 != 0:
        if (hf_pow_ratio > ref_pow_ratio) and (lf_pow_ratio * 100 < 50.0) and (lf_pow_ratio * 100 > 25.0):
            print "NREM 2 detected" # There are 2 peaks one at 12 Hz and another at 15 Hz and the peak is reasonably high (2x the valley)
            stage = 2
    # NREM stages 3 and 4, also known as SWS, is dominated by delta wave (0 - 4 Hz) activity.        
    if (stage == 0) and (lf_pow_ratio * 100 > 40.0):
        stage = 3
        print "NREM 3/4 detected"
    elif(stage == 0):
        print "NREM 1 detected"
        stage = 1
        
    #if(stage != ans):
     #   print "Wrong result, expected " + str(ans)  + " Got " + str(stage)
    
    #return stage
    return stage
    

def find_peak(data):
    
    # search all the samples in the arrany and find the largest peak. Then find its height
    peak = 0;
    peak_i = 0;
#    for i in np.arange(0,len(data)):
#        if data[i] > peak:
#            peak = data[i]
#            peak_i = i     and (data[i-1] > data[i-2])      and (data[i+2] < data[i+1])
    
    for i in np.arange(2,len(data)-2):
        if( (data[i] > data[i-1])  and (data[i+1] < data[i]) ):
            peak = data[i]
            peak_i = i
    
    #print "peak i = " + str(peak_i) #+ " " + str(data)
    valley = peak
    pk_height = 0
    # find the valley before the peak to get its height
    for i in np.arange(0,peak_i):
        #print str(data[peak_i - i-1])
        if data[peak_i - i-1] < valley:
            valley = data[peak_i - i-1]
            pk_height = peak - valley
            #print "valley = " + str(valley)
    if valley != 0:
        pk_height = peak / valley
    print "peak height " + str(pk_height)
    return(pk_height)



def plot_hypnogram(eeg, stages, srate):
    """
    This function takes the eeg, the stages and sampling rate and draws a 
    hypnogram over the spectrogram of the data.
    """
    
    fig,ax1 = plt.subplots()  #Needed for the multiple y-axes
    
    #Use the specgram function to draw the spectrogram as usual
    y_lim = 40;
    plt.specgram(eeg/np.sum(eeg),NFFT=512,Fs=srate)

    #Label your x and y axes and set the y limits for the spectrogram
    ax1.set_ylim((0,y_lim))
    ax1.set_xlim((0,len(eeg)/srate))
    plt.title ('Hypnogram')
    ax1.set_xlabel('Time in Seconds')
    ax1.set_ylabel('Frequency in Hz')
    
    ax2 = ax1.twinx() #Necessary for multiple y-axes
    
    #Use ax2.plot to draw the hypnogram.  Be sure your x values are in seconds
    #HINT: Use drawstyle='steps' to allow step functions in your plot
    ax2.plot(np.arange(0,len(stages))*30,stages,drawstyle='steps')

    #Label your right y-axis and change the text color to match your plot
    ax2.set_ylabel('NREM Stages',color='b')

 
    #Set the limits for the y-axis 
    ax2.set_ylim(0.5,3.5)
    ax2.set_xlim((0,len(eeg)/srate))
    #Only display the possible values for the stages
    ax2.set_yticks(np.arange(1,4))
    
    #Change the left axis tick color to match your plot
    for t1 in ax2.get_yticklabels():
        t1.set_color('b')
    
    #Title your plot    


        
def classifier_tester(classifiedEEG, actualEEG):
    """
    returns percent of 30s epochs correctly classified
    """
    epochs = len(classifiedEEG)
    incorrect = np.nonzero(classifiedEEG-actualEEG)[0]
    percorrect = (epochs - len(incorrect))/epochs*100
    
    print 'EEG Classifier Performance: '
    print '     Correct Epochs = ' + str(epochs-len(incorrect))
    print '     Incorrect Epochs = ' + str(len(incorrect))
    print '     Percent Correct= ' + str(percorrect) 
    print 
    return percorrect
  
    
def test_examples(examples, srate):
    """
    This is one example of how you might write the code to test the provided 
    examples.
    """
    i = 0
    bin_size = 30*srate
    c = np.zeros((4,len(examples[1,:])/bin_size))
    while i + bin_size < len(examples[1,:]):
        for j in range(1,4):
            c[j,i/bin_size] = classify_epoch(examples[j,range(i,i+bin_size)],srate,j)
        i = i + bin_size
    
    totalcorrect = 0
    num_examples = 0
    for j in range(1,4):
        canswers = np.ones(len(c[j,:]))*j
        correct = classifier_tester(c[j,:],canswers)
        totalcorrect = totalcorrect + correct
        num_examples = num_examples + 1
    
    average_percent_correct = totalcorrect/num_examples
    print 'Average Percent Correct= ' + str(average_percent_correct) 
    return average_percent_correct
    

def run_test(eeg,srate):
    i = 0
    bin_size = 30*srate
    stages = np.zeros((len(eeg)/bin_size))
    while i + bin_size < len(eeg):
        stages[i/bin_size] = classify_epoch(eeg[range(i,i+bin_size)],srate)
        i = i + bin_size
    
    
    return stages

def classify_eeg(eeg,srate):
    """
    DO NOT MODIFY THIS FUNCTION
    classify_eeg takes an array of eeg amplitude values and a sampling rate and 
    breaks it into 30s epochs for classification with the classify_epoch function.
    It returns an array of the classified stages.
    """
    bin_size_sec = 30
    bin_size_samp = bin_size_sec*srate
    t = 0
    classified = np.zeros(len(eeg)/bin_size_samp)
    while t + bin_size_samp < len(eeg):
       classified[t/bin_size_samp] = classify_epoch(eeg[range(t,t+bin_size_samp)],srate)
       t = t + bin_size_samp
    return classified
        
##########################
#You can put the code that calls the above functions down here    
if __name__ == "__main__":
    #YOUR CODE HERE
    
    plt.close('all') #Closes old plots.
    
    ##PART 1
    #Load the example data
    examples, srate = load_examples('example_stages.npz')

    x_lim = 30
    #Plot the psds
#    Row 0: REM sleep
#    Row 1: Stage 1 NREM sleep
#    Row 2: Stage 2 NREM sleep
#   Row 3: Stage 3 and 4 NREM sleep (for this problem set, we will combine stages 3 and 4,as in the AASM criteria.)

    print 'Sampling Rate = ' + str(srate)
#    Pxx = np.zeros((4,srate*2+1)) #
#    freq = np.zeros((4,srate*2+1))
#    plt.figure()
#    plt.title('Power Spectral Density')
#    plt.subplot(411)
#    (Pxx[0],freq[0]) = plt.psd(examples[0]/np.sum(examples[0]),NFFT=512,Fs=srate)
#    plt.xlim((0,x_lim))
#    plt.title ('REM')
#    plt.subplot(412)
#    plt.title ('Stage 1 NREM')
#    (Pxx[1],freq[1]) = plt.psd(examples[1]/np.sum(examples[1]),NFFT=512,Fs=srate)
#    plt.xlim((0,x_lim))
#    plt.subplot(413)
#    plt.title ('Stage 2 NREM')
#    (Pxx[2],freq[2]) =plt.psd(examples[2]/np.sum(examples[2]),NFFT=512,Fs=srate)
#    plt.xlim((0,x_lim))
#    plt.subplot(414)
#    plt.title ('Stage 3/4 NREM')
#    (Pxx[3],freq[3]) =plt.psd(examples[3]/np.sum(examples[3]),NFFT=512,Fs=srate)
#    plt.xlim((0,x_lim))
#    plt.show();
    
    
    #Plot the spectrograms
    #plot_example_spectrograms(examples,srate)
    #Test the examples
    #plt.figure()
    #plt.plot(Pxx[2])
    #plt.bar(np.arange(0,srate*2+1),Pxx[2])
    #test_examples(examples,srate)
    #Load the practice data
    #Load the practice answers
    #Classify the practice data
    #Check your performance
    #epoch = 30
    #classify_epoch(examples[3,epoch:epoch+128*30],srate)
    #Generate the hypnogram pls

    eeg, srate = load_eeg('practice_eeg.npz')
    stages = load_stages('practice_answers.npz')
    plot_hypnogram(eeg,stages,srate)
    
    # Hypnogram for test data set
    plt.figure()
    eeg, srate = load_eeg('test_eeg.npz')
    stages = run_test(eeg, srate)
    plot_hypnogram(eeg,stages,srate)
    #plt.specgram(eeg/np.sum(eeg),NFFT=512,Fs=srate)
#    test = classify_eeg(eeg,srate)
#    classifier_tester(test,stages)

    
