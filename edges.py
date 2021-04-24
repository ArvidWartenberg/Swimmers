import numpy as np
from collections import deque
edges = {1:{2,3}, 2:{4}, 5:{6}}

assignment = {}
def recursion(key, cluster):
    if assignment.__contains__(key): return
    assignment[key] = cluster
    if not edges.__contains__(key): return
    connections = edges.pop(key)
    for con in list(connections):
        recursion(con, cluster)

cluster = 0

while bool(edges):
    k = list(edges.keys())[0]
    recursion(k,cluster)
    cluster += 1

clusters = {}
for key in list(assignment):
    if not clusters.__contains__(assignment[key]):
        clusters[assignment[key]] = 1
    else:
        clusters[assignment[key]] += 1

print(assignment)
print(clusters)
print()
