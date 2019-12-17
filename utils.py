import numpy as np
from random import random, seed
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
import datetime as dt
from statistics import mean

def infection_time(event_list, p, seed_node):
    """

    Parameters
    ----------
    event_list
    p
    seed_node

    Returns
    -------
    infection_times
    infection_list

    """

    # seed(12)

    infection_times = {}
    found_seed = False

    for event in event_list:
        if event["Source"] == seed_node and found_seed is False:
            infection_times[event["Source"]] = event["StartTime"]
            found_seed = True

        if event["Source"] in infection_times and event["StartTime"] >= infection_times[
            event["Source"]] and random() <= p:
            if event["Destination"] not in infection_times:
                infection_times[event["Destination"]] = event["EndTime"]
            else:
                if infection_times[event["Destination"]] > event["EndTime"]:
                    infection_times[event["Destination"]] = event["EndTime"]

    infection_list = list(infection_times.values())
    infection_list.sort()
    return infection_times, infection_list


def create_bins(start, end, n_bins):
    """
    Creates a set of linear bins.

    Parameters
    -----------
    start: minimum value of data, int
    end: maximum value of data, int
    n_bins: number of bins, int

    Returns
    --------
    bins: a list of linear bin edges
    """
    bins = np.linspace(start, end, n_bins)
    return bins


def plot_avg_prevalence(infection_times_list, infection_prob, n_nodes, bins):
    """

    Parameters
    ----------
    infection_times_list
    infection_prob
    n_nodes
    bins

    Returns
    -------

    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    dateconv = np.vectorize(dt.datetime.fromtimestamp)
    date = dateconv(bin_centers)
    prevalences = []
    for list, prob in zip(infection_times_list, infection_prob):
        for l in list:
            counts, _, _ = binned_statistic(
                x=l,
                values=l,
                bins=bins,
                statistic='count')

            cum_counts = np.cumsum(counts)
            prevalence = cum_counts / n_nodes
            prevalences.append(prevalence)

        avg_prevalence = np.array(prevalences)
        avg_prevalence = avg_prevalence.mean(0)

        ax.plot(date, avg_prevalence, label=prob)
        prevalences = []
        # ax.get_xaxis().get_major_formatter().set_useOffset(False)
        # ax.get_xaxis().get_major_formatter().set_scientific(False)

    ax.legend()
    plt.suptitle(r'Averaged prevalence of the disease for each of the infection probabilities')
    ax.set_xlabel(r'time')
    ax.set_ylabel(r'averaged prevalence $\rho(t)$')

    fig.autofmt_xdate(bottom=0.2, rotation=20, ha='right')
    fig.savefig("./plots/averaged_prevalence.pdf")
