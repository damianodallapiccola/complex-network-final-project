import numpy as np
from random import random, seed
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic, spearmanr
import datetime as dt
from statistics import mean
import scipy.stats


def infection_time(event_list, p, seed_node, immunized_nodes=None):
    """

    Parameters
    ----------
    immunized_nodes
    event_list
    p
    seed_node

    Returns
    -------
    infection_times
    infection_list

    """

    # seed(12)

    if immunized_nodes is None:
        immunized_nodes = []
    infection_times = {}
    found_seed = False

    for event in event_list:
        if event["Source"] == seed_node and found_seed is False:
            infection_times[str(event["Source"])] = event["StartTime"]
            found_seed = True

        if event["Destination"] not in immunized_nodes and str(event["Source"]) in infection_times and event[
            "StartTime"] >= infection_times[str(event["Source"])] and random() <= p:
            if str(event["Destination"]) not in infection_times:
                infection_times[str(event["Destination"])] = event["EndTime"]
            else:
                if infection_times[str(event["Destination"])] > event["EndTime"]:
                    infection_times[str(event["Destination"])] = event["EndTime"]

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


def plot_avg_prevalence_probs(infection_times_list, infection_prob, n_nodes, bins):
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
    plt.suptitle(r'Averaged prevalence of the disease with different infection probabilities')
    ax.set_xlabel(r'time')
    ax.set_ylabel(r'averaged prevalence $\rho(t)$')

    fig.autofmt_xdate(bottom=0.2, rotation=20, ha='right')
    fig.savefig("./plots/averaged_prevalence_probs.pdf")


def plot_avg_prevalence_nodes(infection_times_list_nodes, seed_nodes_labels, n_nodes, bins):
    """

    Parameters
    ----------
    infection_times_list_nodes
    seed_nodes_labels
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
    for list, label in zip(infection_times_list_nodes, seed_nodes_labels):
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

        ax.plot(date, avg_prevalence, label=label)
        prevalences = []
        # ax.get_xaxis().get_major_formatter().set_useOffset(False)
        # ax.get_xaxis().get_major_formatter().set_scientific(False)

    ax.legend()
    plt.suptitle(r'Average prevalence of the disease with different seed nodes')
    ax.set_xlabel(r'time')
    ax.set_ylabel(r'averaged prevalence $\rho(t)$')

    fig.autofmt_xdate(bottom=0.2, rotation=20, ha='right')
    fig.savefig("./plots/averaged_prevalence_nodes.pdf")


def plot_and_spearman_task4(infection_times_median, clustering_coefficient_net, degree_net, strength_net,
                            betweenness_centrality_net, n_nodes):
    # ordered list of values, the index represent the node
    infection_times_median_list = []
    clustering_coefficient_net_list = []
    degree_net_list = []
    strength_net_list = []
    betweenness_centrality_net_list = []

    for i in range(n_nodes):
        infection_times_median_list.append(infection_times_median[str(i)])
        clustering_coefficient_net_list.append(clustering_coefficient_net[str(i)])
        degree_net_list.append(degree_net[str(i)])
        strength_net_list.append(strength_net[str(i)])
        betweenness_centrality_net_list.append(betweenness_centrality_net[str(i)])

    # dateconv = np.vectorize(dt.datetime.fromtimestamp)
    # date = dateconv(infection_times_median_list)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(clustering_coefficient_net_list, infection_times_median_list, alpha=0.5)
    # plt.ylim(dateconv(dateconv([min(infection_times_median_list)-(min(infection_times_median_list)%(3600*24))])), max(date))
    plt.suptitle(r'Median infection times as a function of the unweighted clustering coefficient')
    ax.set_xlabel(r'clustering coefficient $c$')
    ax.set_ylabel(r'median infection time')
    fig.set_figwidth(6.7)
    fig.savefig("./plots/t4_clustering_coefficient.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(degree_net_list, infection_times_median_list, alpha=0.5)
    plt.suptitle(r'Median infection times as a function of the degree')
    ax.set_xlabel(r'degree $k$')
    ax.set_ylabel(r'median infection time')
    fig.set_figwidth(6.7)
    fig.savefig("./plots/t4_degree_net.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(strength_net_list, infection_times_median_list, alpha=0.5)
    plt.suptitle(r'Median infection times as a function of the strength')
    ax.set_xlabel(r'strength $s$')
    ax.set_ylabel(r'median infection time')
    fig.set_figwidth(6.7)
    fig.savefig("./plots/t4_strength_net.pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(betweenness_centrality_net_list, infection_times_median_list, alpha=0.5)
    plt.suptitle(r'Median infection times as a function of the unweighted betweenness centrality')
    ax.set_xlabel(r'betweenness centrality')
    ax.set_ylabel(r'median infection time')
    fig.set_figwidth(6.7)
    fig.savefig("./plots/t4_betweenness_centrality_net.pdf")

    # Spearman rank-correlation coefficient
    print("Spearman rank-correlation coefficient between median infection time and: ")
    print("- clustering coefficient: " + str(
        spearmanr(infection_times_median_list, clustering_coefficient_net_list).correlation))
    print("- degree: " + str(spearmanr(infection_times_median_list, degree_net_list).correlation))
    print("- strength: " + str(spearmanr(infection_times_median_list, strength_net_list).correlation))
    print("- betweenness centrality: " + str(
        spearmanr(infection_times_median_list, betweenness_centrality_net_list).correlation))


def plot_avg_prevalence_immunization(infection_times_list_immunization, immunization_strategy_labels, n_nodes, bins):
    """

    Parameters
    ----------
    infection_times_list_immunization
    immunization_strategy_labels
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
    for list, label in zip(infection_times_list_immunization, immunization_strategy_labels):
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

        ax.plot(date, avg_prevalence, label=label)
        prevalences = []
        # ax.get_xaxis().get_major_formatter().set_useOffset(False)
        # ax.get_xaxis().get_major_formatter().set_scientific(False)

    ax.legend()
    plt.suptitle(r'Average prevalence of the disease with different immunization strategies')
    ax.set_xlabel(r'time')
    ax.set_ylabel(r'averaged prevalence $\rho(t)$')

    fig.autofmt_xdate(bottom=0.2, rotation=20, ha='right')
    fig.savefig("./plots/t5_averaged_prevalence_immunization.pdf")
