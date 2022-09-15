Overall remarks:
- the osm data preparation is inspired from jafari

-----------------------------------------------------------------------------------------------------------------------

Acquisition
- download the data from OSM

Simplification
- reduce the complexity of intersections, merge nodes within a certain distance

Reconstruct the lane composition
- use tags in osm to generate a string containing the information about individual lanes
- calculate overall stats: number of lanes, total width, oneway y/n

Merge consecutive links
- remove unnecessary nodes

Merge parallel links
- merge links with same u-v combination

Export to SHP for manual corrections

(Do the manual corrections in QGIS)

Reimport the SHP as a graph with geometries

Match with elevation model (pending)
- is this necessary?

Detect dead ends

Match type of intersections (low prio)

Match public transport (pending)
- match ways with route relations
- or by GTFS matching

Match environmental factors (low prio)
- trees, buildings, water, etc. (ask meister what factors are important in his choice model)

Set minimum moto lanes (pending)
- for each way, set which directions must remain (e.g. both for a tramway street or a dead end)
- decide whether the lanes must be separated or one 2-dir lane is fine

Dieting (pending)
- reassign lanes to cycling, ~50% where possible

Build road hierarchy (pending)
- primary rods, secondary, residential, dead-ends, etc.

Split roads into subnetworks (pending)
- use the primary network to split the residential roads into small subnetworks

Reorganize one-way streets (pending)

Generate new lane composition (pending)
- for each way, generate a string showing a possible new lane composition
- calculate overall stats

Generate directed links with capacity for each mode (pending)

Export (pending)
- SHP
- MATSim XML needed?