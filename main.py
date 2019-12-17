import networkx as nx
import pandas as pd
import random
import scipy.stats
from datetime import datetime
import numpy as np
from utils import infection_time, create_bins, plot_avg_prevalence_probs, plot_avg_prevalence_nodes, plot_and_spearman_task4
from si_animator import visualize_si

def main():
    event_data = np.genfromtxt('data/events_US_air_traffic_GMT.txt', names=True, dtype=int)
    event_data.sort(order=['StartTime'])
    network = nx.read_weighted_edgelist('data/aggregated_US_air_traffic_network_undir.edg')
    n_nodes = network.number_of_nodes()

    # creation of bins for the plots
    min_timestemp = min(event_data, key=lambda item:item["StartTime"])[2]
    max_timestemp = max(event_data, key=lambda item:item["EndTime"])[3]
    n_bins = 50
    bins = create_bins(min_timestemp, max_timestemp, n_bins)

    ######################################
    #               task 1               #
    ######################################

    # infection_times, infection_list = infection_time(event_data, 1, 0)
    # print("Node 41 infection time: "+str(infection_times['41'])+" "+str(datetime.fromtimestamp(infection_times['41'])))

    # # animation of the infection
    # visualize_si(np.array(infection_list), save_fname="./simulations/infection_simulation_prob1_seed0.mp4")

    ######################################
    #               task 2               #
    ######################################

    # seed_node = 0
    # infection_prob = [0.01, 0.05, 0.1, 0.5, 1.0]
    # infection_times_list_avg = []
    # infection_times_list_probs = []
    # for prob in infection_prob:
    #     for i in range(10):
    #         _, infection_list = infection_time(event_data, prob, seed_node)
    #         infection_times_list_avg.append(infection_list)
    #     infection_times_list_probs.append(infection_times_list_avg)
    #     infection_times_list_avg = []
    #
    # plot_avg_prevalence_probs(infection_times_list_probs, infection_prob, n_nodes, bins)

    ######################################
    #               task 3               #
    ######################################

    # infection_prob = 0.1
    # seed_nodes = [0, 4, 41, 100, 200]
    # seed_nodes_labels =['ABE', 'ATL', 'ACN', 'HSV', 'DBQ']
    #
    # infection_times_list_avg = []
    # infection_times_list_nodes = []
    # for seed_node in seed_nodes:
    #     for i in range(10):
    #         _, infection_list = infection_time(event_data, infection_prob, seed_node)
    #         infection_times_list_avg.append(infection_list)
    #     infection_times_list_nodes.append(infection_times_list_avg)
    #     infection_times_list_avg = []
    #
    # plot_avg_prevalence_nodes(infection_times_list_nodes, seed_nodes_labels, n_nodes, bins)

    ######################################
    #               task 4               #
    ######################################

    # infection_times_list = []
    # for i in range(50):
    #     infection_prob = 0.5
    #     seed_node = random.randint(0, n_nodes)
    #     infection_times, _ = infection_time(event_data, infection_prob, seed_node)
    #     infection_times_list.append(infection_times)
    #
    # infection_times_df = pd.DataFrame(infection_times_list)
    # infection_times_median = dict(infection_times_df.median())
    # clustering_coefficient_net = nx.clustering(network)
    # degree_net = nx.degree(network)
    # strength_net = nx.degree(network, weight="weight")
    # betweenness_centrality_net = nx.betweenness_centrality(network)
    #
    # plot_and_spearman_task4(infection_times_median, clustering_coefficient_net, degree_net, strength_net, betweenness_centrality_net, n_nodes)






if __name__ == '__main__':
    main()
