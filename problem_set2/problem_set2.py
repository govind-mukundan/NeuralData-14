#
#  NAME
#    problem_set2_solutions.py
#
#  DESCRIPTION
#    Open, view, and analyze action potentials recorded during a behavioral
#    task.  In Problem Set 2, you will write create and test your own code to
#    create tuning curves.
#

#Helper code to import some functions we will use
import numpy as np
import matplotlib.pylab as plt
import matplotlib.mlab as mlab
from scipy import optimize
from scipy import stats
from collections import Counter # count duplicated and stuff


def load_experiment(filename):
    """
    load_experiment takes the file name and reads in the data.  It returns a
    two-dimensional array, with the first column containing the direction of
    motion for the trial, and the second column giving you the time the
    animal began movement during thaht trial.
    """
    data = np.load(filename)[()];
    return np.array(data)

def load_neuraldata(filename):
    """
    load_neuraldata takes the file name and reads in the data for that neuron.
    It returns an arary of spike times.
    """
    data = np.load(filename)[()];
    return np.array(data)
    
def bin_spikes(trials, spk_times, time_bin):
    """
    bin_spikes takes the trials array (with directions and times) and the spk_times
    array with spike times and returns the average firing rate for each of the
    eight directions of motion, as calculated within a time_bin before and after
    the trial time (time_bin should be given in seconds).  For example,
    time_bin = .1 will count the spikes from 100ms before to 100ms after the 
    trial began.
    
    dir_rates should be an 8x2 array with the first column containing the directions
    (in degrees from 0-360) and the second column containing the average firing rate
    for each direction
    
    Note that some trials are repeated. i.e. the same angle may appear more than once
    trials[:,0] = motion angles in radians (the physical activity)
    trials[:,1] = timestamp of each activity
    spk_times[:] = timestamp of each AP recorded for ONE neuron (in Seconds)
    
    dir_rate[:,0] = motion angles in radians
    dir_rate[:,1] = firing rate for that angle in the given window
    """
    angles_dict = Counter(trials[:,0]) # we get a dictionary of the values and their counts
    dir_rates = np.zeros( (len(angles_dict),2 ) )
    angles = angles_dict.items()
    index = 0
    # for each angle sum all the APs over all the trials. angle[0] contains the number of trials for that angle
    for angle in angles: # select a particular angle
        fire_cnt = 0
        for a in range(0,len(trials[:,0])):
            if(angle[0] == trials[a,0]):
                activity_time = trials[a,1]
                for api in range(0,len(spk_times)):
                    if((spk_times[api] >= (activity_time - time_bin)) and (spk_times[api] <= (activity_time + time_bin)) ):
                        fire_cnt = fire_cnt + 1
                        #print "Fire at activity:" + str(activity_time) + "AP Time: " + str(spk_times[api]) + "Angle:" + str(angle[0])
        # Update the (angle, fire count) into the OP array
        # We need to divide by the nunmber of trials to get the average spike count per trial
        # Divide by 2*time_bin to convert the spike count to Firing rate in spikes / second
        dir_rates[index] = [angle[0], fire_cnt /(angle[1]* 2 * time_bin)]
        index = index + 1
    
    dir_rates = dir_rates[dir_rates[:,0].argsort()] # sort by angle
    # argsort() returns the indexes of the sorted elements
    print dir_rates

    # Now lets plot the data
    #plt.figure()
    width = 45
    ax = plt.subplot(2,2,1)
    rects1 = ax.bar(dir_rates[:,0] - width/2, dir_rates[:,1],width)
    ax.set_xlabel("Direction of Motion (degrees)")
    ax.set_ylabel("Firing Rate (spikes/s)")
    ax.set_title("Example Neuron Tuning Curve")
    ax.set_xlim([-width/2,315 + width/2])
    ax.set_xticks(dir_rates[:,0])
        # derive the labels for the x-ticks
    label = []
    for i in range(0,len(dir_rates[:,0])):
        label.append(str(int(dir_rates[i,0])))
        
    ax.set_xticklabels(label)
    
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(11)
    
    # http://matplotlib.org/examples/pylab_examples/polar_demo.html
    # for the Polar plot, duplicate the first value into the value for 360
    #dir_rates = np.append(dir_rates, [360,dir_rates[0,1]])
    theta = np.append(dir_rates[:,0], 360)
    r = np.append(dir_rates[:,1], dir_rates[0,1])
    ax = plt.subplot(222,polar=True)
    ax.set_title("Example Neuron Tuning Curve")
    ax.plot(np.deg2rad(theta),r,label="Firing Rate (spikes/s)")
    ax.legend(loc=8,fontsize=7)

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(11)
        
    plt.show()

    
    
    return dir_rates
    
def plot_tuning_curves(direction_rates, title):
    """
    This function takes the x-values and the y-values  in units of spikes/s 
    (found in the two columns of direction_rates) and plots a histogram and 
    polar representation of the tuning curve. It adds the given title.
    """

    
def roll_axes(direction_rates):
    """
    roll_axes takes the x-values (directions) and y-values (direction_rates)
    and return new x and y values that have been "rolled" to put the maximum
    direction_rate in the center of the curve. The first and last y-value in the
    returned list should be set to be the same. (See problem set directions)
    Hint: Use np.roll()
    """
    # We want to arrange the angles to that the max firing rate is in the middle.
    # So if the list is even, we need to make it odd by adding the first sample into the end
    
    direction_rates = direction_rates[direction_rates[:,0].argsort()] # sort by angle
    
    xs = direction_rates[:,0]
    ys = direction_rates[:,1]    
    id_mid = len(xs) / 2
    count = len(xs)  
    # find the index of the peak firing rate  
    id_maxfr = np.argmax(ys)
    if id_mid >= id_maxfr:
        roll_degrees = id_mid - id_maxfr
    else:
        roll_degrees = id_mid - id_maxfr + count
        
    ys = np.roll(ys,roll_degrees)
    xs = np.roll(xs,roll_degrees)
    
    for i in range(0,roll_degrees):
         xs[i] = xs[i] - 360
    
    # make the first and last sample equal
    if(len(xs) % 2 == 0):
        new_xs = np.append(xs, xs[len(xs) - 1] + 45)
        new_ys = np.append(ys, ys[0] )   
    
    roll_degrees = roll_degrees * 45
    
    print new_xs
    print new_ys
    print roll_degrees
        
    
    return new_xs, new_ys, roll_degrees    
    

def normal_fit(x,mu, sigma, A):
    """
    This creates a normal curve over the values in x with mean mu and
    variance sigma.  It is scaled up to height A.
    """
    n = A*mlab.normpdf(x,mu,sigma)
    return n

def fit_tuning_curve(centered_x,centered_y):
    """
    This takes our rolled curve, generates the guesses for the fit function,
    and runs the fit.  It returns the parameters to generate the curve.
    """
    max_y = np.amax(centered_y)  #What is the biggest y-value? (This estimates the amplitude of the curve,A)
    max_x = centered_x[np.argmax(centered_y)] #Where is the biggest y-value? (This estimates the mean of the curve, mu)
    sigma = 90 #Here we are approximating one standard deviation of our normal distribution (which should be around the width of 2 bars).

    print max_y
    print max_x
    p, cov = optimize.curve_fit(normal_fit,centered_x, centered_y, p0=[max_x, sigma,max_y])

    return p
    


def plot_fits(direction_rates,fit_curve,title):
    """
    This function takes the x-values and the y-values  in units of spikes/s 
    (found in the two columns of direction_rates and fit_curve) and plots the 
    actual values with circles, and the curves as lines in both linear and 
    polar plots.
    """
    ax = plt.subplot(2,2,3)
    ax.hold( True)
    plt.plot(fit_curve[:,0],fit_curve[:,1]) # Plot the smooth curve
    plt.plot(direction_rates[:,0], direction_rates[:,1], 'o') # Plot the actual raw data points
    
    ax.set_xlabel("Direction of Motion (degrees)")
    ax.set_ylabel("Firing Rate (spikes/s)")
    ax.set_title(title)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(11)
    
    ax = plt.subplot(2,2,4,polar=True)
    ax.hold(True)
    theta = np.append(direction_rates[:,0], 360)
    r = np.append(direction_rates[:,1], direction_rates[0,1])
    ax.plot(np.deg2rad(theta),r,'o')
    
    theta = np.append(fit_curve[:,0], 360)
    r = np.append(fit_curve[:,1], fit_curve[0,1])
    ax.plot(np.deg2rad(theta),r,label="Firing Rate (spikes/s)")
    ax.set_title(title)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(11)
    
    ax.legend(loc=8,fontsize=7)

def von_mises_fitfunc(x, A, kappa, l, s):
    """
    This creates a scaled Von Mises distrubition.
    """
    return A*stats.vonmises.pdf(x, kappa, loc=l, scale=s)


    
def preferred_direction(fit_curve):
    """
    The function takes a 2-dimensional array with the x-values of the fit curve
    in the first column and the y-values of the fit curve in the second.  
    It returns the preferred direction of the neuron (in degrees).
    """
    # to find the angle with maxium response, just find the x-value of the max y-value
    id_maxfr = np.argmax(fit_curve[:,1])
    pd = fit_curve[id_maxfr,0]
    print "Preferred Direction of Neuron = " + str(pd) + "Degrees"
    return pd
    

def unroll_axes(x,y,roll_degrees):
    """
    Takes x and y values of the fit curve and unrolls it. It's assumed that the values are spaced by 1
    so roll_degrees of 90 means 90 samples needs to be unrolled
    """
    roll_degrees = roll_degrees
    xs = x[0:len(x) - 1]
    for i in range(0,roll_degrees):
         xs[i] = xs[i] + 360    
    ys = np.roll(y,-roll_degrees)
    xs = np.roll(x,-roll_degrees)
    

    # make the first and last sample equal
    new_xs = np.append(xs, xs[len(xs) - 1] + 1)
    new_ys = np.append(ys, ys[0] ) 
        
    #print new_xs
    return new_xs,new_ys
    
##########################
#You can put the code that calls the above functions down here    
if __name__ == "__main__":
    trials = load_experiment('trials.npy')   
    spk_times = load_neuraldata('neuron2.npy') 
    dir_rate = bin_spikes(trials, spk_times, 0.1)
    cx,cy,rd = roll_axes(dir_rate)
    
    # Now that we have X and Y values arranged in the order of a normal distribution, we can try to fit a curve
    # The idea is that we need to find a curve or function that when given our X values, will return Y values that are close to what we actually got
    p = fit_tuning_curve(cx,cy) # Get the parameters for a normal curve that gives similar results to our Y values
    cx = np.arange(cx[0],cx[-1]) # Get a smooth range of X'es
    fit_ys = normal_fit(cx,p[0],p[1],p[2]) # Find the Y's from the fitted curve
    #plt.figure()
    #plt.plot(cx,fit_ys)
    # Now this curve is rolled, so lets unroll the whole curve so we get back the normal angles
    cx,fit_ys = unroll_axes(cx,fit_ys,rd)

    
    fit_curve = np.empty([len(cx),len(fit_ys)])
    fit_curve[:,0] = cx
    fit_curve[:,1] = fit_ys
    preferred_direction(fit_curve)
    
    plot_fits(dir_rate,fit_curve,"Example Tuning Curve - Fit")
