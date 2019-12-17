import networkx as nx
from datetime import datetime
import numpy as np
from utils import infection_time, create_bins, plot_avg_prevalence


def main():
    event_data = np.genfromtxt('data/events_US_air_traffic_GMT.txt', names=True, dtype=int)
    event_data.sort(order=['StartTime'])
    network = nx.read_weighted_edgelist('data/aggregated_US_air_traffic_network_undir.edg')
    n_nodes = network.number_of_nodes()

    ##########
    # task 1 #
    ##########

    infection_times, infection_list = infection_time(event_data, 1, 0)
    print("Node 41 infection time: " + str(infection_times[41]) +" "+ str(datetime.fromtimestamp(infection_times[41])))

    # animation of the infection
    # visualize_si(np.array(infection_list), save_fname="si_viz_example.mp4")

    ##########
    # task 2 #
    ##########

    infection_prob = [0.01, 0.05, 0.1, 0.5, 1.0]
    infection_times_list_avg = []
    infection_times_list_probs = []
    for prob in infection_prob:
        for i in range(10):
            _, infection_list = infection_time(event_data, prob, 0)
            infection_times_list_avg.append(infection_list)
        infection_times_list_probs.append(infection_times_list_avg)
        infection_times_list_avg = []


    min_timestemp = min(event_data, key=lambda item:item["StartTime"])[2]
    max_timestemp = max(event_data, key=lambda item:item["EndTime"])[3]
    n_bins = 50
    bins = create_bins(min_timestemp, max_timestemp, n_bins)

    plot_avg_prevalence(infection_times_list_probs, infection_prob, n_nodes, bins)






if __name__ == '__main__':
    main()
