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
    (in degrees) and the second column containing the average firing rate
    for each direction
    """
    
    #Set up your bin size
    start_bin = time_bin
    stop_bin = time_bin
    
    #Get just the trial times
    trial_times = trials[:,1]
    
    #Initialize your array for each trial
    spikes_per_trial = np.zeros(len(trial_times))
    
    #Count the number of spikes in the bin around each trial
    for t in np.arange(0,len(trial_times)):
        spikes_per_trial[t] = np.count_nonzero(np.logical_and(spk_times > trial_times[t] - start_bin, spk_times < trial_times[t] + stop_bin))
      

    #Find the directions of motion used in the experiemnt
    directions = np.unique(trials[:,0])
    #Initialize the dir_rates output 
    rates = np.zeros(len(directions))
    for dirs in range(len(directions)):
        d  = directions[dirs]
        indices = plt.find(trials[:,0] == d)  #This lets you group by direction
        rates[dirs] = sum(spikes_per_trial[indices])/len(indices) #average over the trials    
    rates = rates/(start_bin + stop_bin) #Convert to firing rate   
    dir_rates = np.column_stack((directions,rates))
    return dir_rates
    
def plot_tuning_curves(direction_rates, title):
    """
    This function takes the x-values and the y-values  in units of spikes/s 
    (found in the two columns of direction_rates) and plots a histogram and 
    polar representation of the tuning curve. It adds the given title.
    """
    x = direction_rates[:,0]
    y = direction_rates[:,1]
    plt.figure()
    plt.subplot(2,2,1)
    plt.bar(x,y,width=45,align='center')
    plt.xlim(-22.5,337.5)
    plt.xticks(x)
    plt.xlabel('Direction of Motion (degrees)')
    plt.ylabel('Firing Rate (spikes/s)')
    plt.title(title)   
        
        
    
    plt.subplot(2,2,2,polar=True)
    r = np.append(y,y[0])
    theta = np.deg2rad(np.append(x, x[0]))
    plt.polar(theta,r,label='Firing Rate (spikes/s)')
    plt.legend(loc=8)
    plt.title(title)
    
def roll_axes(direction_rates):
    """
    roll_axes takes the x-values (directions) and y-values (direction_rates)
    and return new x and y values that have been "rolled" to put the maximum
    direction_rate in the center of the curve.  It also returns the number of
    degrees that you "rolled" the graph - to use for unrolling late. The first 
    and last y-value in the returned list should be set to be the same. (See problem set directions)
    Hint: Use np.roll()
    """
    x = direction_rates[:,0]
    y = direction_rates[:,1]
    md = np.argmax(y)
    dir_rolled = np.roll(y,4-md)
    roll_degrees = x[4-md]
  
    new_xs = x-roll_degrees 
    new_xs = np.append(new_xs,new_xs[0]+360)
    
    new_ys = np.append(dir_rolled,dir_rolled[0])   
    
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
    max_y = np.amax(centered_y)
    max_x = centered_x[np.argmax(centered_y)]
    sigma = 90

    pn, covn = optimize.curve_fit(normal_fit,centered_x,centered_y,p0=[max_x,sigma,max_y])
    return pn
    
def roll_and_fit(direction_rates):
    """
    Here we run through the steps to roll our axes, fit our curve, and roll
    back to the original axes.  We return an array with one column for the x-values
    and one for the y-values
    """
    new_xs, new_ys, roll_degrees = roll_axes(direction_rates)
    pn = fit_tuning_curve(new_xs,new_ys)
    curve_xs = np.arange(new_xs[0],new_xs[-1])
    fit_ys = normal_fit(curve_xs,pn[0],pn[1],pn[2])  
    final_ys = np.roll(fit_ys,-int(roll_degrees))
    final_xs = np.arange(0,360)    
    return np.column_stack((final_xs,final_ys))

def plot_fits(direction_rates,fit_curve,title):
    """
    This function takes the x-values and the y-values  in units of spikes/s 
    (found in the two columns of direction_rates and fit_curve) and plots the 
    actual values with circles, and the curves as lines in both linear and 
    polar plots.
    """
    plt.subplot(2,2,3)
    plt.plot(direction_rates[:,0],direction_rates[:,1],'o',hold=True)
    plt.plot(fit_curve[:,0],fit_curve[:,1])
    plt.title(title)
    plt.xlabel('Direction of Motion (degrees)')
    plt.ylabel('Firing Rate (spikes/s)')
    plt.xlim(0,360)
    
    plt.subplot(2,2,4,polar=True)
    plt.polar(np.deg2rad(direction_rates[:,0]),direction_rates[:,1],'o')
    plt.polar(np.deg2rad(fit_curve[:,0]),fit_curve[:,1],label='Firing Rate (spikes/s)')
    plt.legend(loc=8)
    plt.title(title)
    

def von_mises_fitfunc(x, A, kappa, l, s):
    """
    This creates a scaled Von Mises distrubition.
    """
    return A*stats.vonmises.pdf(x, kappa, loc=l, scale=s)

def vm_fit(direction_rates):
    """
    This generates the curve fit with a Von Mises distribution.  We return an array  
    with one column for the x-values and one column for the y-values.
    """
    x = direction_rates[:,0]
    y = direction_rates[:,1]
    peak_guess = np.deg2rad(x[np.argmax(y)])
    guesses = [np.amax(y),4.,peak_guess,1]
    p, cov = optimize.curve_fit(von_mises_fitfunc, np.deg2rad(x), y, guesses)
    curve_xs = np.arange(x[0],x[-1]+45)
    fit_ys = von_mises_fitfunc(np.deg2rad(curve_xs),p[0],p[1],p[2],p[3])  
    return np.column_stack((curve_xs,fit_ys)) 
    
def preferred_direction(fit_curve):
    """
    The function takes a 2-dimensional array with the x-values of the fit curve
    in the first column and the y-values of the fit curve in the second.  
    It returns the preferred direction of the neuron (in degrees).
    """
    x = fit_curve[:,0]
    y = fit_curve[:,1]
    pd = x[np.argmax(y)]
    return pd
    
        
##########################

#You can put the code that calls the above functions down here    
if __name__ == "__main__":
    trials = load_experiment('trials.npy')   
    spikes = load_neuraldata('example_spikes.npy') 
   # spikes = load_neuraldata('neuron3.npy')
    direction_rates = bin_spikes(trials,spikes,0.1)
    plot_tuning_curves(direction_rates,'Example Tuning Curve')


    fit_curve = roll_and_fit(direction_rates)
    plot_fits(direction_rates,fit_curve,'Example Tuning Curve - Fit')
    print preferred_direction(fit_curve)

    plot_tuning_curves(direction_rates, 'Example Tuning Curve')
    vmf = vm_fit(direction_rates)
    plot_fits(direction_rates,vmf,'Example Tuning Curve - VM Fit')
    print preferred_direction(vmf)

