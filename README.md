# MDASummer2020
Summer 2020 Internship
This repository consists of the compilation of my work with MD Anderson for Summer 2020
The heatmap files show the step by step data processing and mapping for the microglia specific genes (to run just execute all of the code)--> followed Suerat guided clusternig tutorial
The scale and ordering files have code to: 
-Across all genes take all of the cells that have a read for that gene 
-Scale nonzero values to one another
-Rescale values between 1 and 2 
-Do that for each gene 
-For genes we want to prioritize, make the maximum value 3 and the minimum value 2 (TMEM119)
-For each cell, take each non zero expression and average those
 (that way this is not dependent on reads but on average expression) (sum would prioritize cells with more reads)(averaging all would prioritize cells with more reads) (higher expression average means more prioritized genes expressed or overall high expression) 
-Order cells based on average 
if error encountered: execute line
        SAMPLENAMEmg.data[is.na(SAMPLENAMEmg.data)] <- 0


