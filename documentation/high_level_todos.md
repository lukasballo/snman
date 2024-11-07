In this document, we keep track of long-term, high-level todos and ideas to remember for future major changes.

**Integrate street graph and lane graph closely together:**

- one single api for adding, removing, and editing nodes, edges, and their lanes
- one single data model that allows the ease of editing and compatibility of a street graph 
  but also has the topological consistency and routability of a lane graph

**Understandable edge and node attributes**

- remove unnecessary attributes and name the rest in an understandable way
- introduce scopes, e.g., osm.layer=1

**Generalize the network design**

- introduce rules with semantics, similar to cga in CityEngine
- generalize the inputs that can be provided by the rebuilding regions to the rules used

**Fluid network optimization**

- relax the rigid order of link elimination, allow a more fluid process where different types
  of lanes are added and removed simultaneously
- try out a "random network forest" approach, where the link elimination is done on multiple
  copies of the network in parallel, with some random variation and then the best result is chosen
- placement of elements with priority value, similar to the "badness" value in LaTeX

**Add support for demand data**

- Add a possibility to add an OD matrix for evaluating the total generalized cost

**Create packages for pip and conda**

**Create extensive documentation**

- a more complex case example where all functionalities are used, with images