# ICFP2021

## Rough Strategy
* Divide figure into regions separated by the cut points (and including them). We'll call these "cut components".
* Compute distances and directions ("proximity field") from vertices to the hole boundary, and which ones are outside.
  * Aggregate this information up to the cut components.
* Estimate optimal figure translation and rotation of cut components around cut points.
* Squash and stretch by cut-point-delimited regions.
  * Constraint: The movement on cut points matches between connected cut components.
* Visualize, visualize, visualize.
* Iterate, iterate, iterate.
