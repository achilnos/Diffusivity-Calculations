#PTt script to plot out the P, T, t landscape for experiments. Goal is to genreate a 3D data set where x is T, y is P and z (possibly to be displayed using colored contour graph is the time it takes for the experiment to reach a predetermined thershold)
#for this to work, I need to modify my Parameters file to produce an array of (P,T,vis) and pass to BC_vector. BC_vector must create an array of boundary conditions which is passed to PTt. PTt can iterate through these and use the results to make a new array of (P, T, time to certain threshold). 