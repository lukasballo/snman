Acquisition
- download the data from OSM

Simplification
- reduce intersections into single points

Reconstruct the lane composition
- use tags in osm to reconstruct the lanes of each link
- calculate overall stats: number of lanes, total width, oneway y/n

Merge consecutive links
- remove unnecessary nodes

Merge parallel links
- merge links with same u-v combination

Match with elevation model (pending)

Detect dead ends

Match type of intersections (pending, low prio)

Match public transport
- match ways with route geometries
- detect directions: pt in one or in both ways (pending)
- categorize the pt lines by relevance: tram in tunnel, local bus, night bus, regular on-street...

Match environmental factors (pending, low prio)
- trees, buildings, water, etc. (ask meister what factors are important in his choice model)

Set minimum car lanes (pending)
- for each way, set which directions must remain (e.g. both for a tramway street or a dead end)
- decide whether the lanes must be separated or one 2-dir lane is fine

Dieting (pending)
- reassign lanes to cycling

Build road hierarchy
- primary rods, secondary, residential, dead-ends, etc.

Split roads into subnetworks (pending)
- use the primary network to split the residential roads into small subnetworks

Reorganize one-way streets (pending)
- One-Way Traffic Organization Problem
- Goal: Creating a reasonable network, similar to how a planner would do it
- Optimization problem, the cost function could be based on a detour measure and/or betweenness centrality

Generate new lane composition (pending)
- for each way, generate a string showing a possible new lane composition
- calculate overall stats

Generate directed links with capacity for each mode + environmental factors for cycling (pending)

Export (pending)
- MATSim XML