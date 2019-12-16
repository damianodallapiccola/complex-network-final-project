import networkx as nx
import matplotlib as mpl
import numpy as np
import scipy
from random import random
from si_animator import visualize_si

def infection_time(event_list, p, seed_node):
    """
    inputs:
    - network
    - event list
        format: Source, Destination, StartTime, EndTime, Duration
    - probability p
    - the seed node id

    outputs:
    - the number of infected nodes as function of time
    - a dictionary of the infection times of with nodes as keys.
    """

    infection_times = {}
    found_seed = False

    for event in event_list:
        if event["Source"] == seed_node and found_seed is False:
            infection_times[event["Source"]] = event["StartTime"]
            found_seed = True

        if event["Source"] in infection_times and event["StartTime"] >= infection_times[event["Source"]] and random() <= p:
            if event["Destination"] not in infection_times:
                infection_times[event["Destination"]] = event["EndTime"]
            else:
                if infection_times[event["Destination"]] > event["EndTime"]:
                    infection_times[event["Destination"]] = event["EndTime"]

    infection_list = list(infection_times.values())
    infection_list.sort()
    return infection_times, infection_list


def main():
    event_data = np.genfromtxt('events_US_air_traffic_GMT.txt', names=True, dtype=int)
    event_data.sort(order=['StartTime'])

    ##########
    # task 1 #
    ##########

    infection_times, infection_list = infection_time(event_data, 1, 0)
    print("Node 41 infection time: " + str(infection_times[41]))

    # animation of the infection
    visualize_si(np.array(infection_list), save_fname="si_viz_example.mp4")

    ##########
    # task 2 #
    ##########





    # network = nx.read_weighted_edgelist('aggregated_US_air_traffic_network_undir.edg')


if __name__ == '__main__':
    main()
