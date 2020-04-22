"""MGAB Custom Time Series Generator

This script allows the user to to create his/her own Mackey Glass 
Anomaly Benchmark. The main function of this module is 
generate_benchmark(args). All parameters are passed to this function 
through a Python dictionary. It is possible to pass an empty dictionary 
or no argument at all. Typically, one would specify a subset of the 
required parameters in the dictionary; the function would then use 
the default values for the remaining parameters.

Usage:
    import mgab
    mgab = generate_benchmark(args)

This file can also be imported as a module and contains the following
functions:

    * print_error - Writes error messages to stderr
    * get_default_params - Generates a dictionary with a set of 
      default parameters
    * generate_benchmark - Generates a MGAB according to the specifications 
      of the user. A list of MG time series with a certain number of 
      anomalies is created. The created time series can be directly written 
      to CSV-files and/or returned by this function and processed further.
    * create_anomalous_time_series - Randomly selects different locations in
      a given time series where anomalies will be inserted
    * set_anomalies - Based on a given list of locations, this function will
      insert anomalies by removing segments of the time series.
    * lyapunov_MLE - Estimates the Maximum Lyapunov Exponent (MLE) of a given 
      time series. A positive MLE is usually considered as an indictation for
      chaotic behaviour in the time series
    * generate_mackey_glass - Generates a Mackey-Glass (MG) time series according
      to the requirements which are passed to the function in a dictionary. 
      This function will run a DDE solver (jitcdde) to obtain the MG time series.
      
The main data structure which is passed to many of the functions is a dictionary,
where all possible options can be specified:
args: Dictionary, containing all user-specified parameters for generating a MGAB 
benchmark. It is not necessary to pass any argument to generate_benchmark(). In 
this case, the original MGAB data is created which was also uploaded to this 
repository. The function will use the pre-computed MG time series from the file 
./data/mgts_len=5000000tau=18n=10.0bet=0.25gam=0.1h=0.9T=1.npy to generate the 
10 original time series. It is also possible to create the whole benchmark from 
scratch, by adjusting the dictionary entry of 'reproduce_original_mgab', 
as described below. 

Currently, the following options are supported:
    - 'verbosity' (int): Integer, describing the level of verbosity. 
         - 0: No outputs to stdout 
         - 1: Print standard Info-messages
         - Larger than 1: Debug Messages 
         + default:  *1*
    - 'output_dir' (str): Directory, in which the resulting files containing the individual 
         time series will be written. The option None allows to turn off the saving. Per default, 
         the files are named 1.csv, 2.csv, ... and have the same format as the original files 
         described above. 
         + default:  './mgab/'
    - 'output_force_override' (bool): If a file with the same name already exists in the output 
         directory, then this will *not* be overwritten per default. This can be changed by setting 
         this option to *True*. So with this option you can force the generator to override benchmark 
         files in the output directory, if a file with that name already exists.
         + default:  *False*
    - 'output_format'(str): Currently, only 'csv' is supported. All output files, saved as 1.csv, 
         2.csv, ..., are comma-seperated-value (CSV) files.
         + default:  *'csv'*
    - 'num_series'(int): The number of individual anomalous time series, which should be generated. 
         Creating many time series might take a lot of time, since the DDE solver has to pre-compute 
         (and if no already pre-computed MG time series is available, see 'mg_ts_path_load') one long 
         MG time series which is then split into the specified number.   
         + default:  *10*
    - 'series_length'(int): The length of one anomalous time series. Similarly to 'num_series', it might 
         be computationally expensive, if too large values are chosen.
         + default:  *100000*
    - 'num_anomalies'(int) : The number of anomalies which are placed into each time series. E.g., if 10 
         time series are generated and this option is set to 10, then the overall benchmark will contain 
         100 anomalies.
         + default:  *10*,
    - 'noise'(str) :   In order to further increase the complexity of the benchmark and to be more similar 
         to real-world problems, there is a possibility to add random noise to the individual anomalous MG 
         time series. It is also possible to turn off the noise component. The amount of noise is controlled
         by 'noise_param'. There are 4 options: 
         + 'rnd_uniform': Adds random uniform noise to the time series.
         + 'rnd_walk': Samples from a random uniform distribution and computes a random walk (cumulative sum) 
            over these values.
         + 'rnd_normal': Sample from a Gaussian normal distribution.
         + None: Do not add any noise to the time series
         + default:  *'rnd_uniform'*
    - 'noise_param'(tuple,float) : For 'rnd_uniform' and 'rnd_walk', specifies the lower and upper bound for 
         the random uniform distribution, respectively. For 'rnd_normal', the first element describes the mean 
         (loc) and the second element the standard deviation (scale).
         + default:  *(-0.01, 0.01)* for 'rnd_uniform', *(-0.001, 0.001)* for 'rnd_walk' 
           and *(0, 0.01)* for 'rnd_normal'.
    - 'min_anomaly_distance'(int) : The minimum distance between 2 anomalies. The difficulty of the benchmark 
         usually increases, if this distance is reduced.
         + default:  *2000*,
    - 'mg_ts_path_load'(str): Full path to a pre-computed MG time series, which can be used to produce the anomalous 
         time series. The file must contain a numpy array ('.npy'-file). Usually, a lot of time can be saved, if a 
         pre-computed MG time series is available, since the DDE solver does not have to be called in this case. However, 
         this time series has to be long enough, to be split into 'num_series' new series of length 'series_length'. 
         The length should be roughly 1.5 x 'num_series' x 'series_length', since also segments of the time series are 
         removed in order to create the anomalies. If this pre-computed time series is too short, an exception will be 
         thrown. In this case, one could reduce 'num_series' or 'series_length' or set this parameter to *None*. If this 
         option is set to *None*, no pre-computed MG time series will be loaded from disk. Instead, the DDE solver will 
         generate a time series which is suffiently long. 
         + default:  *None* 
    - 'mg_tau'(float) :  Parameter *τ* in the MG equation.
         + default:  *18.0*
    - 'mg_n'(float):  Parameter *n* in the MG equation.
         + default:  *10.0*
    - 'mg_beta'(float) : Parameter *β* in the MG equation.  
         + default:  *0.25*
    - 'mg_gamma'(float) : Parameter *γ* in the MG equation.
         + default:  *0.1*
    - 'mg_history'(float) : Value for the constant history *h*, which is required by the DDE solver as initial condition.
         + default:  *0.9*
    - 'mg_T'(float): Step-size parameter T for the DDE solver. Usually, 'mg_T=1' is sufficient. Smaller step-sizes will 
         reduce the "gaps" between the evaluated points of the MG equation and increase the number of data points which 
         are generated.
         + default:  *1.0*
    - 'mg_ts_dir_save'(str): In order to possibly save time in the future, this option can be used to specify a file 
         (only the directory) where the pre-computed MG time series shall be saved. The filename will contain all necessary 
         parameters, which allow the user to later find a certain setting. Since the filename is unique (according to the 
         setting), no check is performed to ensure that no duplicate file is present in the folder. If a file with the same 
         filename should already be present, then this file should contain exactly the same data as the generated one.    
         + default:  *None*,
    - 'seed'(int): In order to allow to reproduce certain settings, the user may specify a seed. If no seed is provided 
         (*None*) then the algorithm will take the value i for the i-th time series as seed (0,1,...). So, if a completly 
         random result is necessary, the user should set a sufficiently random seed (e.g., a timestamp in ms, etc.)
         + default: *None* 
    - 'min_length_cut'(int): This option specifies the minimum size of the segments which are removed from the time series in 
         order to create an anomaly. It does not make sense to make trivial cuts (e.g. cut segment of length 1, which might 
         happen if adjacent points have the largest similarity in a certain range). Usually, the default value is a good choice.
         + default:  *100*,
    - 'max_sgmt'(int):  Maximum segment length in which we want to find the most similar values. In order to create an anomaly 
         we compare the values (and the derivatives) of the time series in two windows with each other. Then, the segment between
         the 2 most similar values is removed and the 2 remaining ends are "stiched" together again: 
         # xxxxxxxxxxx [window1] xxxxx [window2] xxxxxxxxxxx. 'max_sgmt' basically describes the size of window1 and window2 
         (which both have the same size).  
         + default:  *100* 
    - 'anomaly_window'(int): Window size of anomaly window which is put around each anomaly. Smaller windows increase the difficulty, 
         since algorithms have to locate the anomalies more accurately.
         + default:  *400*
    - 'order_derivative'(int): Until which (numerical) derivative do we want to compare the similarity of the points? (0->only value, 
         1-> value and 1st derivative, ...)
         + default:  *3*
    - 'reproduce_original_mgab': This option is used to generate the original MGAB, which is described in the beginning of this page 
         and for which the data has been added to this repository. If this option is set to a value which is not *None*, then most of 
         the previous options will be ignored and the default settings will be taken, so that the original results can be re-produced. 
         Currently, it is only possible to adjust the other options 'output_dir', 'output_force_override' and 'output_format' (which 
         currently also only has one possibility).  There are 3 options:
         + 'generate_new_mg': With this option, the whole benchmark is re-computed from scratch. This can take a long time, since the 
            DDE solver has to be run again.
         + 'use_precomputed_mg': Using this option, the pre-computed MG time series in 
            ./data/mgts_len=5000000tau=18n=10.0bet=0.25gam=0.1h=0.9T=1.npy will be used to generate the benchmark. The compuation time 
            should be limited in this case. 
         + None: This value is usually chosen, since in most cases we do not want to re-produce the old time series, but rather create our own new benchmark.
         + default:  *None*
"""


def print_error(err_msg):
    """Prints an error message to stderr 

        Parameters
        ----------
        err_msg : str
            The error message, which is written to stderr.

        Raises
        ------
    """
    import sys
    sys.stdout.flush()  # first get rid of all other messages
    sys.stderr.write(err_msg)
    sys.stderr.write("\n")
    sys.stderr.flush() # Make sure, the error is immediately shown
# --------------------------------------------------------------------------------------------------------------

def get_default_params():
    """ Returns the default options which can be used to generate a new 
        MGAB with the standard settings.

        Parameters
        ----------
        None
        
        Raises
        ------
 
        Returns
        -------
        a dictionary containing all default options which are required to generate
        the MGAB. Note that these options are used in the generation process, if
        not explicitely specified by the user. The options of the dictionary are
        described above.
    """
    return {
        'verbosity': 1,  # 0: No outputs, 1: Standard Info-messages, >=2: Debug Messages
        'output_dir': 'mgab',
        'output_force_override': False,
        # With this parameter you can force the generator to override benchmark files in the output directory, if a file with that name already exists
        'output_format': 'csv',  # 'csv' (currently, only csv supported)
        'num_series': 10,
        'series_length': 100000,
        'num_anomalies': 10,
        'noise': 'rnd_uniform',  # 'rnd_normal', 'rnd_walk', None
        'noise_param': (-0.01, 0.01),  # range for 
        'min_anomaly_distance': 2000,
        'mg_ts_path_load': None,
        # Default path to a long pre-computed MG time series, which can be used to generate the time series.

        # Parameters for the generation of the Mackey-Glass Time Series:
        'mg_tau': 18.0,
        'mg_n': 10.0,
        'mg_beta': 0.25,
        'mg_gamma': 0.1,
        'mg_history': 0.9,
        'mg_T': 1,
        'mg_ts_dir_save': None,

        # Details for anomly insertion process:
        'seed': None,  # None: Take the value i for the i-th time series as seed (0,1,...)
        'min_length_cut': 100,  # It does not make sense to make a trivial cut (e.g. cut segment of length 1)
        'max_sgmt': 100,  # maximum segment length in which we want to find the closest value
        'anomaly_window': 400,
        # window size of anomaly window which we put around the anomaly (TODO: divide by 2 in function)
        'order_derivative': 3,
        # until which derivative do we want to compare the similarity of the points (0->only value, 1-> value and 1st derivative, ...)
        'reproduce_original_mgab': None  # None, 'generate_new_mg', 'use_precomputed_mg'
    }
# --------------------------------------------------------------------------------------------------------------

def generate_benchmark(args: dict = {'reproduce_original_mgab': 'use_precomputed_mg'}):
    """ Generates a MGAB according to the specifications of the user. A list of MG time series with a 
    certain number of anomalies is created. The created time series can be directly written to 
    CSV-files and/or returned by this function and processed further.

        Parameters
        ----------
        args : dict, optional
            A dictionary where different options can be set by the user. Note that it is not required
            to set any options. In this case the defaults will be used. Also note that it makes a
            difference if an empty dictionary is passed or if no argument is passed at all to this
            function. If an empty dictionary is passed, a benchmark similar to the original MGAB 
            will be created, however the generated time series will not be exactly the same. In order
            to exactly reproduce the time series of the original MGAB just call the function without
            any argument, e.g.:
                original_mgab = mgab.generate_benchmark()
        
        Raises
        ------
            None
            However, several exceptions might be thrown (and caught) within this function which will 
            lead to error messages. We try to handle most errors to ensure that no results are lost
            after a possibly lenghty generation process. Ideally, even if some files cannot be saved,
            the function should return a list with the generated time series when the generation 
            procedure is completed.
            Furthermore, the function aborts if inconsistent options are found. Since these checks
            are done in the beginning, no computation time should be lost.
 
        Returns
        -------
        A list of pandas DataFrames. Each DataFrame contains one anomalous time series and also the 
        anomaly labels. In total there are 4 columns: index, value, is_anomaly, is_ignored. The last
        column can be used to allow a warm-up phase for algorithms and not to count wrong detections
        in the inital 256 time steps.
    """
    import sys, os
    import numpy as np
    import time

    par = get_default_params()  # Load default params first

    # Check, if there are some unrecognized arguments
    for k in args.keys():
        if k not in par:
            print_error("Error: Unrecognized argument: " + str(k) + ". Stop here!")
            return None
    par.update(args)  # overwrite the defaults with user-provided params

    # If we want to reproduce the original benchmark, then reload the defaults again
    if par['reproduce_original_mgab'] is not None:
        par = get_default_params()
        # Only allow a few user-specified parameters
        par['reproduce_original_mgab'] = args['reproduce_original_mgab']
        if 'output_dir' in args: par['output_dir'] = args['output_dir']
        if 'output_force_override' in args: par['output_force_override'] = args['output_force_override']
        if 'output_format' in args: par['output_format'] = args['output_format']
        # Furthermore, set the time series path to the one in the directory data, if wanted
        if par['reproduce_original_mgab'] == 'generate_new_mg':
            par['mg_ts_path_load'] = None
        else:
            par['mg_ts_path_load'] = "./data/mgts_len=5000000tau=18n=10.0bet=0.25gam=0.1h=0.9T=1.npy"

    # Check, if the specified lengths make sense
    # How long does the time series have to be, if we have this many anomalies in this setting?
    min_ts_length = par['num_anomalies'] * (
            par['min_anomaly_distance'] * 4) // 3  # x1.33 should be sufficient in most cases
    if par['series_length'] < min_ts_length:
        print_error("The specified length of the series is too short for the current setup!")
        print_error("Consider increasing 'series_length' or reducing 'num_anomalies' or 'min_anomaly_distance'. ")
        print_error("Stopping now!")
        return None

    # if noise_params is not specified and the noise should be Gaussian, then adjust loc and scale
    if ('noise_params' not in args) and par['noise'] == 'rnd_normal':
        par['noise_params'] = (0, 0.01)  # loc=0 and scale=0.01

    verbose = par['verbosity']
    if verbose > 0:
        print("Generating the MGAB with the following parameters: ")
        print(par)
        print(
            "Note that some parameters might have been set to the default values, if these were not specified by you!")
        print()

    if (par['output_dir'] is None):
        if verbose > 0:
            print("No output directory specified. The benchmark files will not be saved.")
            print("However, generate_benchmark() will still return a list with the anomalous time series")
    else:
        # If output directory does not exist, try to create it here:
        if not os.path.isdir(str(par['output_dir'])):
            if verbose > 0:
                print("Directory does not yet exist. I will try to create it!")
            try:
                # Create target Directory
                os.makedirs(str(par['output_dir']))
                if verbose > 0:
                    print("Directory ", str(par['output_dir']), " created ")
            except:
                print_error(
                    "Error: Something went wrong. Could not create the output folder " + str(par['output_dir']) +
                    "\n. I will still continue with the next steps and return the resulting time series!\n")
                import traceback
                traceback.print_exc()
                par[
                    'output_dir'] = None  # Make sure, we dont try to write anyting into the non-existant directory later

    # First, we have to either load a pre-computed MG time series or we have to generate it here
    mg_required_length = par["num_series"] * par["series_length"]  # This is the minimum length
    # Since we remove segments later, we have to consider this as well: 
    # In the extreme case we will remove a range from X until X + min_length_cut + 2*max_sgmt
    mg_required_length += par["num_series"] * par["num_anomalies"] * (
            par['min_length_cut'] + 2 * par['max_sgmt'])  # par['min_anomaly_distance'] + 
    mg_required_length = (mg_required_length * 6) // 5  # just to make sure
    if par['reproduce_original_mgab'] is not None:
        mg_required_length = 15 * 10**4 * par["num_series"]  # This was used for the original benchmark
        if verbose > 0:
            print(
                "Rebuilding original benchmark. For that I need a MG time series of overall length:",
                mg_required_length)
    # round up the length of each time series to next full 10-thousand in this strange way. This is done to be able to reproduce
    # the original benchmark
    # mg_required_length = int(np.ceil(mg_required_length/10000*1.0/par["num_series"]))*10000*par["num_series"] 

    if par['mg_ts_path_load'] is not None:
        if verbose > 0:
            print("Loading a pre-computed MG time series: ", par['mg_ts_path_load'])
            print("Note: Loading a pre-computed MG time series will ignore the parameters mg_* !")
            print(
                "If you require other settings for your MG time series, you have to generate a new MG time series and remove the 'mg_ts_path_load' parameter!")
            print()
        long_mg_series = np.load(par['mg_ts_path_load'],
                                 allow_pickle=False)  # Will throw an exception, if file was not found
    else:
        # Compute the MG time series here
        # The required length depends on a few parameters 
        gen_args = par.copy()
        gen_args["series_length"] = mg_required_length  # Overwrite the length parameter, since we need some longer TS
        if verbose > 0:
            rough_time_prediction = round((8.8 / 483600.0) * mg_required_length * 1.25, 1)
            print("Generating a new MG time series of length", mg_required_length, "using jitcdde.")
            print("Note that this might take a long time.")
            print("For this case we estimate a computation time of at least", rough_time_prediction, "minutes!")
            print("...")

        start = time.time()
        long_mg_series = generate_mackey_glass(gen_args)
        end = time.time()
        if verbose > 0:
            print("Finished Generating the MG time series! Wall-clock time required (in minutes):",
                  round((end - start) / 60.0, 1))
            print()
        # If we want to save the generated time series, then we have to check, if the directory exsists
        if par['mg_ts_dir_save'] is not None:
            if os.path.isdir(str(par['mg_ts_dir_save'])):
                file_name = "/mgts_" + "len=" + str(mg_required_length) + "tau=" + str(par['mg_tau']) + "n=" + str(
                    par['mg_n']) + "bet=" + str(args['mg_beta']) + "gam=" + str(par['mg_gamma']) + "h=" + str(
                    par['mg_history']) + "T=" + str(par['mg_T']) + ".npy"
                full_path = par['mg_ts_dir_save'] + file_name
                if verbose > 0:
                    print("Saving MG time series to the file:", full_path)
                try:
                    np.save(full_path, long_mg_series, allow_pickle=False)
                except:
                    print_error(
                        "Error: Something went wrong. Could not save the generated MG time series. Will still continue with the next steps...\n")
                    import traceback
                    traceback.print_exc()
            else:
                print_error("Error: Could not find the specified directory " + str(
                    par['mg_ts_dir_save']) + "\n.Will still continue with the next steps...\n")
    if long_mg_series.shape[0] < mg_required_length:
        print_error("Error: The length of the pre-computed MG time series was too short. Required: " + str(
            mg_required_length) + ". Actual:" + str(long_mg_series.shape[0]) + "\nStop here!\n")
        return None
    else:
        long_mg_series = long_mg_series[:mg_required_length]  # this is done to be able to reproduce the original MGAB

    # Split the long time series into chunks of equal size
    length = long_mg_series.shape[0] // par["num_series"]
    if verbose > 2:
        print("long_mg_series.shape[0]:", long_mg_series.shape[0], "par['num_series']:", par["num_series"])
    mg_list = []
    for i in range(par["num_series"]):
        mg_list.append(long_mg_series[i * length: (i + 1) * length])

    # Now we can start creating the individual time series of the benchmark
    all_ts = []
    np.random.seed(par["seed"])  # Set the seed, if an external one was given!
    for i in range(par["num_series"]):
        seed = par["seed"]
        if seed is None:  # Use time series number as seed
            seed = i
            np.random.seed(seed)
        ts_i = create_anomalous_time_series(mg_list[i], par)
        all_ts.append(ts_i)

        # save the benchmark time series into the specified directory
        if par['output_dir'] is not None:
            file_path = par['output_dir'] + "/" + str(i + 1)
            if par['output_format'] == '':  # No other type supported yet!
                try:
                    if os.path.isfile(file_path + ".npy") and not par['output_force_override']:
                        print_error("A file with the name " + str(file_path + ".npy") + " already exists.")
                        print_error(
                            "Consider removing the file or allowing to override the file with the parameter 'output_force_override'!")
                        raise Exception("File already exists!")
                    np.save(file_path + ".npy", ts_i, allow_pickle=False)
                except:
                    import traceback
                    traceback.print_exc()
                    print_error(
                        "Something went wrong while saving the file " + file_path + ".npy! I will not try to save any more files.\nBut I will still attempt to continue and return the results!.")
                    par['output_dir'] = None  # prevent that more files are saved
            else:
                if par['output_format'] != 'csv':
                    print_error("Unknown File Format specified. Will just save as CSV!\n")
                try:
                    if os.path.isfile(file_path + ".csv") and not par['output_force_override']:
                        print_error("A file with the name " + str(file_path + ".csv") + " already exists.")
                        print_error(
                            "Consider removing the file or allowing to override the file with the parameter 'output_force_override'!")
                        raise Exception("File already exists!")
                    ts_i.to_csv(file_path + ".csv")
                except:
                    import traceback
                    traceback.print_exc()
                    print_error(
                        "Something went wrong while saving the file " + file_path + ".csv! I will not try to save any more files.\nBut I will still attempt to continue and return the results.")
                    par['output_dir'] = None  # prevent that more files are saved

    # return a list of the individual time series (as pandas DataFrames)
    np.random.seed(None)  # reset the seed
    return all_ts
# --------------------------------------------------------------------------------------------------------------

def create_anomalous_time_series(series, args: dict):
    """ Generates a MGAB according to the specifications of the user. A list of MG time series with a 
    certain number of anomalies is created. The created time series can be directly written to 
    CSV-files and/or returned by this function and processed further.

        Parameters
        ----------
        series : np.array
            A one-dimensional array, containing a (normal) time series in which anomalies will be 
            inserted. This array has to be longer than the final length which is specified in the
            args dictionary. This is due to the fact that the anomaly insertion process removes
            several segments of the time series.
        args : dict
            A dictionary containing all required options. Refer to the module description above to 
            see a detailed discussion of all options. Not all options will be required in this function
            but it also does not harm to pass non-used options.
        
        Raises
        ------
        Exception
            If the selection of the anomaly locations fails (since it loops until a valid selection
            is found), then a general Exception is thrown. This Exception might be thrown if too many
            anomalies have to be inserted in a too short time series. In this case there might be no
            way to insert all anomalies and also ensuring a minimum distance between the anomalies.
            
 
        Returns
        -------
        A time series with the specified number of anomalies as pandas DataFrame. If requested, 
        also a certain amount of noise is added to the resulting time series.
        
    """
    import numpy as np
    verbose = args["verbosity"]
    length = args["series_length"]
    num_anomalies = args["num_anomalies"]
    min_anomaly_distance = args["min_anomaly_distance"]

    if verbose > 1:
        import matplotlib.pyplot as plt
        print("series.shape:", series.shape)
        plt.figure(figsize=(20, 6))
        tmp = np.random.randint(length, size=None)
        plt.plot(series[tmp:tmp + 1000])
        plt.show()

    loop_counter = 0
    while True:
        anomaly_positions = np.random.randint((length - min_anomaly_distance * num_anomalies) / num_anomalies,
                                              size=num_anomalies) + min_anomaly_distance
        anomaly_positions = int(0.95 * length) - np.cumsum(anomaly_positions)
        if np.all(anomaly_positions > 0):
            break
        loop_counter += 1
        if loop_counter > 100:
            raise Exception("Something went wrong generating the anomalies!")
    anomaly_positions.sort()

    new_series, plots = set_anomalies(series=series, anomalies_idx=anomaly_positions, args=args)
    if verbose > 1:
        print("create_anomalous_time_series(): new_series.shape: ", new_series.shape)

    noise = 0.0
    if (args['noise'] == 'rnd_uniform') or (args['noise'] == 'rnd_walk'):
        noise = np.random.uniform(low=args["noise_param"][0], high=args["noise_param"][1],
                                  size=new_series["value"].shape)  # random uniform
        if args['noise'] == 'rnd_walk':
            noise = noise.cumsum()
    elif (args['noise'] == 'rnd_normal'):
        noise = numpy.random.normal(loc=args["noise_param"][0], scale=args["noise_param"][1],
                                    size=new_series["value"].shape)
    if verbose > 2:
        print("new_series.shape:", new_series.shape, "noise.shape:", noise.shape)
    new_series["value"] += noise

    return new_series.iloc[0:length]
# --------------------------------------------------------------------------------------------------------------

def lyapunov_MLE(args: dict, length=10000):
    """Calculates the Maximum Lyapunov Exponent (MLE) for a given setting 

        If the argument length isn't passed in, the default length of
        10k is used.

        Parameters
        ----------
        length : int, optional
            The length of the MG time series which is generated to estimate the MLE
        args : dict
            A dictionary with all necessary parameters describing the MG equation
            The following options have to be specified in the dictionary:
               * "mg_tau"
               * "mg_n"
               * "mg_beta"
               * "mg_gamma"
               * "mg_history"
               * "step_size"

        Raises
        ------
        
        Returns
        -------
        The estimated Maximum Lyapunov Exponent (MLE). Positive values are usually
        taken as an indication that 
    """

    from jitcdde import jitcdde_lyap, y, t
    import numpy
    from scipy.stats import sem

    τ = args["mg_tau"]
    n = args["mg_n"]
    β = args["mg_beta"]
    γ = args["mg_gamma"]
    history = args["mg_history"]
    step_size = args["step_size"]

    f = [β * y(0, t - τ) / (1 + y(0, t - τ) **  n) - γ * y(0)]
    n_lyap = 1
    DDE = jitcdde_lyap(f, n_lyap=n_lyap)

    DDE.set_integration_parameters(atol=1.0e-16, rtol=1.0e-5, min_step=1.0e-10)
    DDE.constant_past([history])
    DDE.step_on_discontinuities()
    data = []
    lyaps = []
    weights = []
    for time in numpy.arange(DDE.t, DDE.t + length, step_size):
        state, lyap, weight = DDE.integrate(time)
        data.append(state)
        lyaps.append(lyap)
        weights.append(weight)

    lyaps = numpy.vstack(lyaps)
    Lyaps = list()
    stderrs = list()
    for i in range(n_lyap):  # we only have one here
        Lyap = numpy.average(lyaps[400:, i], weights=weights[400:])
        stderr = sem(lyaps[400:, i])  # Note that this is only an estimate
        print("%i. Lyapunov exponent: % .4f +/- %.4f" % (i + 1, Lyap, stderr))
        Lyaps.append(Lyap)
        stderrs.append(stderr)
    return (Lyaps[0], stderrs[0])  # return the MLE (Maximum Lyapunov Exponent)
# --------------------------------------------------------------------------------------------------------------


def generate_mackey_glass(args: dict):
    """ Generates a new Mackey-Glass time series using an accurate DDE solver (jitcdde).
        Depending on the length of the series, the solver can require quite some time.
        We observed that for a series of length 100k about 3 minutes are needed. However,
        this value is only a rough hint, since other parameter settings might require more
        time and also different target platforms / operating systems have large impact on
        the computation time.
    
        Parameters
        ----------
        args : dict
           The following options have to be specified in the dictionary:
               * "series_length"
               * "mg_tau"
               * "mg_n"
               * "mg_beta"
               * "mg_gamma"
               * "mg_history"
               * "step_size"
            
        
        Raises
        ------
       
        Returns
        -------
        The generated Mackey-Glass time series as a numpy array of shape (series length,)
        
        
    """
    from jitcdde import jitcdde, y, t
    import numpy as np
    τ = args["mg_tau"]
    n = args["mg_n"]
    β = args["mg_beta"]
    γ = args["mg_gamma"]
    history = args["mg_history"]
    stepsize = args["mg_T"]  # fixed for the moment (check, if smaller values are required)
    length = args["series_length"]

    f = [β * y(0, t - τ) / (1 + y(0, t - τ) **  n) - γ * y(0)]
    DDE = jitcdde(f)
    DDE.set_integration_parameters(atol=1.0e-16, rtol=1.0e-10)  # min_step = 1.0e-15

    DDE.constant_past([history])
    # DDE.step_on_discontinuities()
    DDE.integrate_blindly(0.0)  # This gives the results comparable to the MATLAB dde23 solver

    data = []
    for time in np.arange(DDE.t, DDE.t + length, stepsize):
        data.append(DDE.integrate(time))

    return np.array(data).squeeze()
# --------------------------------------------------------------------------------------------------------------


def set_anomalies(series, anomalies_idx: list, args: dict):
    """ Inserts anomalies into a time series by removing segments at certain positions. 
        If the split points, where the segments will be removed, are chosen carefully,
        the manipulation will be later hardly visible to the human eye. We use an approach,
        where we compare the values and derivatives in two windows with each other. Those
        two points (one from each window) which have the largest similarity will be selected
        as split points and the segment in between will be removed. Then the two remaining
        ends of the time series are "stiched" together again.
    
        Parameters
        ----------
        series: np.array
           A numpy array containing the normal (nominal) time series in which the anomalies 
           are to be inserted
        anomalies_idx: list
           A list of locations, where the anomalies will be inserted. Note that these locations
           are only rough positions, since we will look in the neighborhood of the anomaly location 
           for suitable points which are very similar. So the actual anomaly might be more on the
           left or the right of the given location.
        args : dict 
           A dictionary containing the necessary options for the insertion process. It does not harm
           to pass unrelated options as well, since they will be ignored. But the following options
           must be present:
           * 'min_length_cut'
           * 'max_sgmt'
           * 'anomaly_window'
           * 'order_derivative'
           * 'verbosity'
           A detailed description of all options can be found in the module description at the beginning
           of this file.
        Raises
        ------
       
        Returns
        -------
        Returns a time series, where the anomalies have been inserted according to the
        specified locations.
        
        
    """
    import numpy as np
    import pandas as pd
    min_length_cut = args[
        'min_length_cut']  # It does not make sense to make a trivial cut (e.g. cut segment of length 1)
    max_sgmt = args['max_sgmt']  # maximum segment length in which we want to find the closest value
    anomaly_window = args['anomaly_window'] // 2  # window size of anomaly window which we put around the anomaly
    order_derivative = args[
        'order_derivative']  # until which derivative do we want to compare the similarity of the points (0->only value, 1-> value and 1st derivative, ...)
    verbose = args['verbosity']

    if verbose > 1:  # Since we will generate some plots
        import matplotlib.pyplot as plt

    real_anomalies_idx = list()
    anomalies_idx = np.sort(anomalies_idx)

    plots = []
    for ad_idx in anomalies_idx:
        gradients = list()
        gradients.append(series)
        for i in range(order_derivative):
            gradients.append(np.gradient(gradients[-1], axis=-1))
        if len(series.shape) > 1:
            all_grads = np.hstack(gradients).T
        else:
            all_grads = np.stack(gradients)

        if verbose > 2:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(20, 6))
            plt.plot(all_grads[0][ad_idx - 50:ad_idx + 50], label='mg1')
            plt.plot(all_grads[1][ad_idx - 50:ad_idx + 50], label='d/dx mg1')
            plt.plot(all_grads[2][ad_idx - 50:ad_idx + 50], label='d^2 / dx^2 mg1')
            plt.legend()
            plt.show()

        mod_window_mean = all_grads.mean(axis=-1, keepdims=True)
        mod_window_std = all_grads.std(axis=-1, keepdims=True)

        # Compare window1 to window2:
        # xxxxxxxxxxx [window1] xxxxx [window2] xxxxxxxxxxx
        mod_window1 = all_grads[:, ad_idx:ad_idx + max_sgmt]
        mod_window2 = all_grads[:, ad_idx + max_sgmt + min_length_cut: ad_idx + min_length_cut + 2 * max_sgmt]

        similarities = np.sqrt(
            (np.apply_along_axis(lambda x: mod_window1 - x[:, np.newaxis], axis=0, arr=mod_window2) **  2).sum(axis=0))

        best_points = np.argwhere(similarities == np.min(similarities))[0]  # (index first window, index 2nd window)
        if verbose > 3:
            print("best_points:", best_points)
            print("ad_idx:", ad_idx)
        idx_first = best_points[0] + ad_idx
        idx_second = best_points[1] + ad_idx + max_sgmt + min_length_cut

        # Cut out a segment
        new_series = np.concatenate([series[:idx_first], series[idx_second:]], axis=0)
        idx_anomaly = idx_first
        real_anomalies_idx.append(idx_anomaly)
        if verbose > 1:
            import matplotlib.pyplot as plt
            print("First cut:", idx_first)
            print("Second cut:", idx_second)
            plt.figure(figsize=(10, 6))
            plt.plot(series, label='original')
            plt.plot(new_series, label='manipulated')
            plt.axvline(x=idx_first, color='r', alpha=0.9, linewidth=.8)
            plt.legend(loc="best")
            plt.xlabel("time / s")
            plt.ylabel("y(t)")
            plt.xlim((idx_first - 50, idx_first + 300))
            plt.show()

        plots.append({"original": series.copy(), "manipulated": new_series.copy(), "idx_cut_first": idx_first,
                      "idx_cut_second": idx_second})
        series = new_series

    if verbose > 1:
        print("All anomalies:", real_anomalies_idx)

    col = ["value"]
    if len(series.shape) > 1:
        col = ["value" + str(i + 1) for i in range(series.shape[-1])]
    series = pd.DataFrame(series, columns=col)
    series["is_anomaly"] = 0
    series["is_ignored"] = 0

    # Ignore anomalies in the very beginning
    series.loc[0:256, ("is_ignored")] = 1

    # set anomaly windows
    for i in real_anomalies_idx:
        series.loc[i - anomaly_window:i + anomaly_window, ("is_anomaly")] = 1

    return series, plots
# --------------------------------------------------------------------------------------------------------------
