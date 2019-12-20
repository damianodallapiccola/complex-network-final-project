import networkx as nx
import pandas as pd
import random
import scipy.stats
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
from utils import infection_time, create_bins, plot_avg_prevalence_probs, plot_avg_prevalence_nodes, plot_and_spearman_task4, plot_avg_prevalence_immunization, infection_edges, plot_scatterplot
from si_animator import visualize_si, plot_network_usa
from collections import Counter

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

    # # ----- task 4 and 5 ----- #
    # clustering_coefficient_net = nx.clustering(network)
    # degree_net = nx.degree(network)
    # strength_net = nx.degree(network, weight="weight")
    # betweenness_centrality_net = nx.betweenness_centrality(network)
    # # ------------------------ #

    # infection_prob = 0.5
    # infection_times_list = []
    # for i in range(50):
    #     seed_node = random.randint(0, n_nodes)
    #     infection_times, _ = infection_time(event_data, infection_prob, seed_node)
    #     infection_times_list.append(infection_times)
    #
    # infection_times_df = pd.DataFrame(infection_times_list)
    # infection_times_median = dict(infection_times_df.median())
    #
    # plot_and_spearman_task4(infection_times_median, clustering_coefficient_net, degree_net, strength_net, betweenness_centrality_net, n_nodes)

    ######################################
    #               task 5               #
    ######################################

    # # nodes immunized
    # imm_neighbour = []
    # range_nodes = set(range(0, n_nodes))
    # while len(imm_neighbour) < 10:
    #     rand_node = random.choice(list(range_nodes))
    #     rand_neighbour = random.choice(list(network.neighbors(str(rand_node))))
    #     if(int(rand_neighbour) not in imm_neighbour):
    #         imm_neighbour.append(int(rand_neighbour))
    #
    #
    # imm_random_node = []
    # range_nodes = set(range(0, n_nodes))
    # for i in range(10):
    #     rand_node = random.choice(list(range_nodes))
    #     imm_random_node.append(rand_node)
    #     range_nodes.remove(rand_node)
    #
    #
    # imm_clustering_coefficient = []
    # d = Counter(clustering_coefficient_net)
    # for k, _ in d.most_common(10):
    #     imm_clustering_coefficient.append(int(k))
    #
    # imm_degree = []
    # highest_degree = sorted(degree_net, key=lambda x: x[1], reverse=True)[:10]
    # for k, _ in highest_degree:
    #     imm_degree.append(int(k))
    #
    # imm_strength = []
    # highest_strength = sorted(strength_net, key=lambda x: x[1], reverse=True)[:10]
    # for k, _ in highest_strength:
    #     imm_strength.append(int(k))
    #
    # imm_betweenness_centrality = []
    # d = Counter(betweenness_centrality_net)
    # for k, _ in d.most_common(10):
    #     imm_betweenness_centrality.append(int(k))
    #
    #
    # # create a set of all the immunized nodes
    # imm_nodes = set(imm_neighbour) | set(imm_random_node) | set(imm_clustering_coefficient) | set(imm_degree) | set(imm_strength) | set(imm_betweenness_centrality)
    # range_seed = set(range(0, n_nodes)) - imm_nodes
    #
    # # extract the seed nodes from a set of nodes not part of the immunized ones
    # seed_nodes = []
    #
    # for i in range(20):
    #     rand_seed = random.choice(list(range_seed))
    #     seed_nodes.append(rand_seed)
    #     range_seed.remove(rand_seed)
    #
    #
    # immunized_nodes_list = []
    # immunized_nodes_list.append(imm_neighbour)
    # immunized_nodes_list.append(imm_random_node)
    # immunized_nodes_list.append(imm_clustering_coefficient)
    # immunized_nodes_list.append(imm_degree)
    # immunized_nodes_list.append(imm_strength)
    # immunized_nodes_list.append(imm_betweenness_centrality)
    # immunization_strategy_labels =['random neighbour', 'random node', 'clustering coefficient', 'degree', 'strength', 'betweenness centrality']
    # infection_prob = 0.5
    #
    #
    # infection_times_list_avg = []
    # infection_times_list_immunization = []
    # for immunized_nodes, imm_strategy in zip(immunized_nodes_list, immunization_strategy_labels):
    #     print(imm_strategy)
    #     for seed_node in seed_nodes:
    #         _, infection_list = infection_time(event_data, infection_prob, seed_node, immunized_nodes)
    #         infection_times_list_avg.append(infection_list)
    #     infection_times_list_immunization.append(infection_times_list_avg)
    #     infection_times_list_avg = []
    #
    # plot_avg_prevalence_immunization(infection_times_list_immunization, immunization_strategy_labels, n_nodes, bins)



    ######################################
    #               task 6               #
    ######################################

    # id_data = np.genfromtxt('data/US_airport_id_info.csv', delimiter=',', dtype=None, names=True)
    # xycoords = {}
    # for row in id_data:
    #     xycoords[str(row['id'])] = (row['xcoordviz'], row['ycoordviz'])
    #
    # edge_list = []
    # for edge in network.edges():
    #     if int(edge[0]) > int(edge[1]):
    #         edge = (edge[1], edge[0])
    #     edge_list.append(edge) # edge_list created to maintain the right order
    #
    # infection_prob = 0.5
    # infecting_edges_fraction = []
    # for i in range(20):
    #     seed_node = random.randint(0, n_nodes)
    #     infecting_edges = infection_edges(event_data, infection_prob, seed_node, edge_list)
    #     infecting_edges_fraction.append(infecting_edges)
    #
    # # calculation of the fraction of times that each link is used for infecting the disease from the results of 20 runs
    # infecting_edges_fraction = (np.sum(np.array(infecting_edges_fraction), 0)/20).tolist()
    #
    # # print Transmission links - fraction
    # fig, ax = plot_network_usa(network, xycoords, edges=edge_list, linewidths=infecting_edges_fraction)
    #
    # plt.suptitle(r'Transmission links ($f_{ij}$)')
    #
    # fig.savefig("./plots/t6_map_fraction.pdf")
    #
    # # print Transmission links - mst
    # maximum_spanning_tree = nx.maximum_spanning_tree(network)
    # fig, ax = plot_network_usa(maximum_spanning_tree, xycoords, edges=list(maximum_spanning_tree.edges))
    #
    # plt.suptitle(r'Transmission links (maximal spanning tree)')
    #
    # fig.savefig("./plots/t6_map_mst.pdf")
    #
    # link_weights = nx.get_edge_attributes(network, 'weight')
    # link_betweenness_centrality = nx.edge_betweenness_centrality(network)
    #
    # # ordered lists (following the order of edge_list)
    # link_weights_list = []
    # link_betweenness_centrality_list = []
    # for edge in edge_list:
    #     if edge in link_weights:
    #         link_weights_list.append(link_weights[edge])
    #     else:
    #         link_weights_list.append(link_weights[(edge[1], edge[0])])
    #     if edge in link_betweenness_centrality:
    #         link_betweenness_centrality_list.append(link_betweenness_centrality[edge])
    #     else:
    #         link_betweenness_centrality_list.append(link_betweenness_centrality[(edge[1], edge[0])])
    #
    #
    # # scatter plot of the transmission fraction as a function of the link weight
    # fig, ax = plot_scatterplot(link_weights_list, infecting_edges_fraction)
    # plt.suptitle(r'Transmission fraction as a function of the link weight')
    # ax.set_xlabel(r'link weight $w_{ij}$')
    # ax.set_ylabel(r'transmission fraction $f_{ij}$')
    # fig.savefig("./plots/t6_scatter_weight.pdf")
    #
    # # scatter plot of the transmission fraction as a function of the link betweenness centrality
    # fig, ax = plot_scatterplot(link_betweenness_centrality_list, infecting_edges_fraction)
    # plt.suptitle(r'Transmission fraction as a function of the link betweenness centrality')
    # ax.set_xlabel(r'unweighted link betweenness centrality $eb_{ij}$')
    # ax.set_ylabel(r'transmission fraction $f_{ij}$')
    # fig.savefig("./plots/t6_scatter_bet_centr.pdf")
    #
    #
    # # Spearman rank-correlation coefficient
    # print("Spearman rank-correlation coefficient between transmission fraction and: ")
    # print("- link weight: " + str(
    #     spearmanr(link_weights_list, infecting_edges_fraction).correlation))
    # print("- betweenness centrality: " + str(
    #     spearmanr(link_betweenness_centrality_list, infecting_edges_fraction).correlation))



if __name__ == '__main__':
    main()
