'''This module provides basic cluster processing for seismic events.'''

from __future__ import print_function
import collections
import numpy as num

epsilon = 1e-6
km = 1000.


class DBEvent:
    '''
    Event condition with respect to the DBSCAN clustering.

    :param event_id: event index (sequential number)
    :param event_type: type of the event according to clustering
    :param directly_reachable: list of ids of events directly reachable
    '''

    def __init__(self, event_id, event_type, directly_reachables):
        self.event_id = event_id
        self.event_type = event_type
        self.directly_reachables = directly_reachables


class ClusterError(Exception):
    pass


def dbscan(simmat, nmin, eps, ncluster_limit):
    '''
    Apply DBSCAN algorithm, reading a similarity matrix and returning a list of
    events clusters

    :param simmat: similarity matrix (numpy matrix)
    :param nmin: minimum number of neighbours to define a cluster
    :param eps: maximum distance to search for neighbours
    '''
    CORE = 1
    NOT_CORE = 2
    EDGE = 3
    ISLE = 4
    REACHED = 5

    nev = len(simmat)
    clusterevents = []
    for i in range(nev):
        event_id = i
        directly_reachables = []
        for j in range(nev):
            if simmat[i, j] <= eps:
                if (i != j):
                    directly_reachables.append(j)

        if len(directly_reachables) >= nmin:
            event_type = CORE
        else:
            event_type = NOT_CORE

        clusterevents.append(
            DBEvent(event_id, event_type, directly_reachables))

    for i in range(nev):
        for j in range(nev):
            if (i in clusterevents[j].directly_reachables) \
                    and (clusterevents[j].event_type == CORE):

                if clusterevents[i].event_type == NOT_CORE:
                    clusterevents[i].event_type = EDGE

        if clusterevents[i].event_type == NOT_CORE:
            clusterevents[i].event_type = ISLE

    eventsclusters = num.zeros((nev), dtype=int)
    actualcluster = -1
    for i in range(nev):
        if clusterevents[i].event_type == ISLE:
            eventsclusters[i] = -1
        else:
            if clusterevents[i].event_type == CORE:
                actualcluster = actualcluster + 1
                reachedevents = []
                pointingevents = []

                eventsclusters[i] = actualcluster
                reachedevents.append(i)
                pointingevents.append(i)

                while len(pointingevents) > 0:
                    newpointingevents = []
                    for j in pointingevents:
                        for k in clusterevents[j].directly_reachables:
                            if clusterevents[k].event_type == CORE:
                                reachedevents.append(k)
                                newpointingevents.append(k)
                                eventsclusters[k] = actualcluster
                                clusterevents[k].event_type = REACHED
                            elif clusterevents[k].event_type == EDGE:
                                reachedevents.append(k)
                                eventsclusters[k] = actualcluster
                                clusterevents[k].event_type = REACHED

                    pointingevents = []
                    for newpointingevent in newpointingevents:
                        pointingevents.append(newpointingevent)

    n_clusters = actualcluster + 1
#   resorting data clusters (noise remains as -1)

    clustersizes = []
    for icl in range(n_clusters):
        lencl = len([ev for ev in eventsclusters if ev == icl])
        clustersizes.append([icl, lencl])

    resorted_clustersizes = sorted(
        clustersizes, key=lambda tup: tup[1], reverse=True)

    resorting_dict = {}
    resorting_dict[-1] = -1

    for icl in range(n_clusters):
        if ncluster_limit is None or icl < ncluster_limit:
            resorting_dict[resorted_clustersizes[icl][0]] = icl
        else:
            resorting_dict[resorted_clustersizes[icl][0]] = -1

    for ievcl, evcl in enumerate(eventsclusters):
        eventsclusters[ievcl] = resorting_dict[evcl]

    return eventsclusters


def get_clusters(events, eventsclusters):
    '''
    Provided a list of events and a list of cluster labels of equal size
    returns a dictionary of clusters with all their events

    i.e.

         cluster -1 -> [event1, event4, eventx]
         cluster  0 -> [event2, event3, eventn]
    '''

    origclusters = {}
    if len(events) != len(eventsclusters):
        raise ClusterError(
            'number of events different from number of cluster labels: %s, %s'
            % len(events), len(eventsclusters))

    for iev, evcl in enumerate(eventsclusters):
        if evcl not in origclusters:
            origclusters[evcl] = [events[iev]]
        else:
            origclusters[evcl].append(events[iev])

    return collections.OrderedDict(sorted(origclusters.items()))
