# BigWhamLink
## a Wolfram LibraryLink interface to BigWham using LTemplate

We ship LTemplate for convenience. 

# Instructions
  - You need to edit the following file: BigWhamLink/BigWhamLink/Kernel/BuildSettings.m and modify the mathematica variable $BigWhamDirectory 
to reflect the path to the BigWham folder (on my machine "/Users/bricelecampion/Documents/Work/Geomechanics/Codes/BigWham/")
  - You may need to edit the path to the mkl and tbb libraries also in BuildSettings.m 
  - Upon first load via << BigWhamLink`, compilation of a layered c++ interface is performed, note that the BigWham static libray must have been built before.

# Examples
- Look at BigWhamLink.nb to get started
- additional examples can be found under the Examples folder
