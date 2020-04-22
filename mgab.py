# ............................................................................
# TODO. Documentation
# Currently refer to https://github.com/MarkusThill/MGAB/blob/master/README.md
# .............................................................................
"""MGAB Time Series Generator

This script allows the user to to create his/her own Mackey Glass 
Anomaly Benchmark

This file can also be imported as a module and contains the following
functions:

    * print_error - 
    * get_default_params
    * generate_benchmark
    * create_anomalous_time_series
    * lyapunov_MLE
    * generate_mackey_glass
"""


def print_error(err_msg):
    import sys
    sys.stdout.flush()  # first get rid of all other messages
    sys.stderr.write(err_msg)
    sys.stderr.write("\n")
    sys.stderr.flush()


# Do not change the default params. If you need other default parameter, you can write a new function... 
def get_default_params():
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
        # Default path to a long pre-computed MG time series, which can be used to generate the time series. Should be a numpy array.

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


def generate_benchmark(args: dict = {'reproduce_original_mgab': 'use_precomputed_mg'}):
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
        mg_required_length = 15 * 10 ** 4 * par["num_series"]  # This was used for the original benchmark
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


def create_anomalous_time_series(series, args: dict):
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

        If the argument `length` isn't passed in, the default length of
        10k is used.

        Parameters
        ----------
        length : int, optional
            The length of the MG time series which is generated to estimate the MLE
        args : dict
            A dictionary with all necessary parameters describing the MG equation

        Raises
        ------
        NotImplementedError
            Not yet.
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

    f = [β * y(0, t - τ) / (1 + y(0, t - τ) ** n) - γ * y(0)]
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
    from jitcdde import jitcdde, y, t
    import numpy as np
    τ = args["mg_tau"]
    n = args["mg_n"]
    β = args["mg_beta"]
    γ = args["mg_gamma"]
    history = args["mg_history"]
    stepsize = args["mg_T"]  # fixed for the moment (check, if smaller values are required)
    length = args["series_length"]

    f = [β * y(0, t - τ) / (1 + y(0, t - τ) ** n) - γ * y(0)]
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
            (np.apply_along_axis(lambda x: mod_window1 - x[:, np.newaxis], axis=0, arr=mod_window2) ** 2).sum(axis=0))

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
