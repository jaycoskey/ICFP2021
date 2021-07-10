# ICFP2021

## Rough Strategy
* Divide figure into regions separated by the cut points (and including them). We'll call these "cut components".
* Call a cut component "free" if it contains only one cut point. It's free to rotate about that cut point.
* Iterate approximate of best (a) rotation of free cut components, (b) figure translation, (c) squash & stretch,
  * Use proximity info:
    * Vertices:
      * Compute distances and directions ("proximity field") from vertices to the hole boundary,
        and which ones are outside.
      * Possibly aggregate this info up to the cut components.
    * Edges: 
      * Compute info on proximity of figure edges to the hole boundary, and find crossings.
* Note: Manipulate the graph of cut components (connected at the cut points),
  noting that some rotations around cut points will impact multiple cut components.
* Visualize, visualize, visualize.
* Iterate, iterate, iterate.
