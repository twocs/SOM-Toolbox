% SOM Toolbox 2.1
% 
% Copyright 1997-2000 by
% Esa Alhoniemi, Johan Himberg, Juha Parhankangas and Juha Vesanto
% Contributed files may contain copyrights of their own.
% 
% SOM Toolbox comes with ABSOLUTELY NO WARRANTY; for details
% see License.txt in the program package. This is free software,
% and you are welcome to redistribute it under certain conditions;
% see License.txt for details.
% 
% 
% Demos
% 
%       demo/som_demo1   SOM Toolbox demo 1: basic properties
%       demo/som_demo2   SOM Toolbox demo 2: basic usage
%       demo/som_demo3   SOM Toolbox demo 3: visualization
%       demo/som_demo4   SOM Toolbox demo 4: data analysis
%       demo/gtm_demo1   GTM functionality demo 1: basic properties
%       demo/gtm_demo2   GTM functionality demo 2: basic usage
% 
% Creation of structs
% 
%          som/som_set   create & set (& check) values to structs
%         som/som_info   print out information on a given struct  
%  som/som_data_struct   create & initialize a data struct 
%   som/som_map_struct   create & initialize a map struct 
% som/som_topol_struct   create & initialize a topology struct 
% som/som_train_struct   create & initialize a train struct 
%     som/som_clstruct   create a cluster struct
%        som/som_clset   set properties in a cluster struct
%        som/som_clget   get stuff from a cluster struct
% 
% Struct conversion and file I/O
% 
%       som/som_vs1to2   converts a version 1.0 struct to version 2.0 struct
%       som/som_vs2to1   converts a version 2.0 struct to version 1.0 struct
%    som/som_read_data   reads a (SOM_PAK format) ASCII data file
%   som/som_write_data   writes a SOM_PAK format codebook file
%    som/som_write_cod   writes a SOM_PAK format data file
%     som/som_read_cod   reads a SOM_PAK format codebook file
% 
% Data preprocessing
% 
%    som/som_normalize   normalize data set
%  som/som_denormalize   denormalize data set 
% som/som_norm_variable   (de)normalize one variable
%       som/preprocess   preprocessing GUI
% 
% Initialization and training functions
% 
%         som/som_make   create, initialize and train a SOM
%     som/som_randinit   random initialization algorithm
%      som/som_lininit   linear initialization algorithm
%     som/som_seqtrain   sequential training algorithm
%   som/som_batchtrain   batch training algorithm
%          som/som_gui   SOM initialization and training GUI
%   som/som_prototrain   a simple version of sequential training: easy to modify
% 
% Clustering algorithms
% 
%       som/som_kmeans   k-means algorithm (was earlier kmeans)
%  som/kmeans_clusters   try and evaluate several k-means clusterings
%       som/neural_gas   neural gas vector quantization algorithm
%      som/som_linkage   hierarchical clustering algorithms
%    som/som_cllinkage   hierarchical clustering of SOM
%   som/som_dmatminima   local minima from distance (or U-) matrix
% som/som_dmatclusters   distance (or U-) matrix based clustering
%     som/som_clspread   spreads clusters to unassinged map units
%       som/som_cldist   calculate distances between clusters
%     som/som_gapindex   gap validity index of clustering
%         som/db_index   Davies-Bouldin validity index of clustering  
% 
% Supervised/classification algorithms
% 
%   som/som_supervised   supervised SOM algorithm
%             som/lvq1   LVQ1 algorithm
%             som/lvq3   LVQ3 algorithm
%              som/knn   k-NN classification algorithm 
%          som/knn_old   k-NN classification algorithm (old version)
% 
% SOM error measures
% 
%      som/som_quality   quantization and topographic error of SOM
%   som/som_distortion   SOM distortion measure
%  som/som_distortion3   elements of the SOM distortion measure
% 
% Auxiliary functions
% 
%         som/som_bmus   calculates BMUs for given data vectors
%     som/som_eucdist2   pairwise squared euclidian distances between vectors
%        som/som_mdist   calculates pairwise distances between vectors 
%       som/som_divide   extract subsets of data based on map
%        som/som_label   give labels to map units
%    som/som_label2num   rcodes string data labels to interger class labels 
%    som/som_autolabel   automatically labels the SOM based on given data
%  som/som_unit_coords   calculates coordinates in output space for map units
%   som/som_unit_dists   distances in output space between map units
%  som/som_unit_neighs   units in 1-neighborhood for each map unit
% som/som_neighborhood   calculates neighborhood matrix for the given map
%    som/som_neighbors   calculates different kinds of neighborhoods 
%       som/som_neighf   calculates neighborhood function values
%       som/som_select   GUI for manual selection of map units
% som/som_estimate_gmm   create Gaussian mixture model on top of SOM
%  som_probability_gmm   evaluate Gaussian mixture model
%      som/som_ind2sub   from linear index to subscript index 
%      som/som_sub2ind   from subscript index to linear index
%      som/som_ind2cod   from linear index to SOM_PAK linear index 
%      som/som_cod2ind   from SOM_linear index to SOM_PAK linear index 
%         som/nanstats   mean, std and median which ignore NaNs
% som/som_modify_dataset   add, remove, or extract samples and components
%     som/som_fillnans   fill NaNs in a data set based on given SOM
%        som/som_stats   statistics of a data set
%       som/som_drmake   calculate descriptive rules for a cluster
%       som/som_dreval   evaluate descriptive rules for a cluster
%     som/som_drsignif   rule significance measures
% 
% Using SOM_PAK from Matlab
% 
%  som/som_sompaktrain   uses SOM_PAK to train a map
%       som/sompak_gui   GUI for using SOM_PAK from Matlab
%      som/sompak_init   call SOM_PAK's initialization programs from Matlab
%  som/sompak_init_gui   GUI for using SOM_PAK's initialization from Matlab
% som/sompak_rb_control   an auxiliary function for sompak_*_gui functions.
%    som/sompak_sammon   call SOM_PAK's Sammon program from Matlab
% som/sompak_sammon_gui   GUI for using SOM_PAK's Sammon program from Matlab
%     som/sompak_train   call SOM_PAK's training program from Matlab
% som/sompak_train_gui   GUI for using SOM_PAK's training program from Matlab 

% Visualization
% 
%         som/som_show   basic visualization
%     som/som_show_add   add labels, hits and trajectories
%   som/som_show_clear   remove extra markers
%   som/som_recolorbar   refresh/reconfigure colorbars
%     som/som_show_gui   GUI for using som_show and associated functions
%         som/som_grid   visualization of SOM grid
%       som/som_cplane   component planes and U-matrices
%     som/som_barplane   bar chart visualization of map
%     som/som_pieplane   pie chart visualization of map
%    som/som_plotplane   plot chart visualization of map
%   som/som_trajectory   launches a GUI for presenting comet-trajectories 
%   som/som_dendrogram   visualization of clustering tree
%   som/som_plotmatrix   pairwise scatter plots and histograms
% som/som_order_cplanes   order and visualize the component planes
%       som/som_clplot   plots of clusters (based on cluster struct)
% som/som_projections_plot   projections plots (see som_projections)
%   som/som_stats_plot   plots of statistics (see som_stats)
% 
% Auxiliary functions for visualization
% 
%             som/hits   calculates hits, or sum of values for each map unit
%         som/som_hits   calculates the response of data on the map
%         som/som_umat   calculates the U-matrix
%              som/cca   curvilinear component analysis projection algorithm
%          som/pcaproj   principal component projection algorithm
%           som/sammon   Sammon's mapping projection algorithm
%   som/som_connection   connection matrix for map 
%   som/som_vis_coords   map unit coordinates used in visualizations
%    som/som_colorcode   create color coding for map/2D data
%     som/som_bmucolor   colors of the BMUs from a given map color code
%    som/som_normcolor   simulate indexed colormap
% som/som_clustercolor   color coding which depends on clustering structure
%  som/som_kmeanscolor   color coding according to k-means clustering
% som/som_kmeanscolor2   a newer version of the som_kmeanscolor function
%   som/som_fuzzycolor   a fuzzy color coding 
%     som/som_coloring   a SOM-based color coding 
%  som/som_projections   calculates a default set of projections
% 
% Report generation stuff
% 
% som/som_table_struct   creates a table struct
% som/som_table_modify   modifies a table struct
%  som/som_table_print   print a table in various formats
%        som/rep_utils   various utilities for printing report elements
%  som/som_stats_table   a table of data set statistics
% som/som_stats_report   report on data set statistics
% 
% Low level routines used by visualization functions
% 
%        som/vis_patch   defines hexagonal and rectangular patches
% som/vis_som_show_data   returns UserData and subplot handles stored by som_show.m
%    som/vis_valuetype   used for type checks 
%     som/vis_footnote   adds a movable text to the current figure 
%      som/vis_trajgui   the actual GUI started by som_trajectory.m 
% som/vis_PlaneAxisProperties   set axis properties in visualization functions
% som/vis_footnoteButtonDownFcn   callback function for vis_footnote.m
% som/vis_planeGetArgs   converts topol struct to lattice, msize argument pair
% som/vis_show_gui_comp   internal function used by som_show_gui.m
% som/vis_show_gui_tool   internal function used by som_show_gui.m 
% 
% DIJKSTRA
% 
% dijkstra/dijkstraKay.m   Shortest paths from nodes 's' to nodes 't' using Dijkstra algorithm
%     dijkstra/LICENCE   MIT License
% 
% GTM
% 
%        gtm/consist.m
%          gtm/dist2.m
%         gtm/eigdec.m
%            gtm/gmm.m
%       gtm/gmmactiv.m
%        gtm/gmmpost.m
%        gtm/gmmprob.m
%            gtm/gtm.m
%          gtm/gtmem.m
%         gtm/gtmem2.m
%       gtm/gtmemseq.m
%   gtm/gtmexpimpute.m
%        gtm/gtminit.m
%       gtm/gtminit2.m
%       gtm/gtmlmean.m
%       gtm/gtmlmode.m
%         gtm/gtmmag.m
%   gtm/gtmmapimpute.m
%        gtm/gtmpost.m
%        gtm/gtmprob.m
%    gtm/gtmvariance.m
%       gtm/gtm_make.m
%       gtm/gtm_show.m
%            gtm/pca.m
%            gtm/rbf.m
%         gtm/rbffwd.m
%       gtm/rbfjacob.m
%       gtm/rbfprior.m
%       gtm/rbfsetfw.m
%       gtm/rbfunpak.m
%
% Contributed information
% 
%       contrib/README
%
% Contributed package GMLVQ
%
% contrib/gmlvq/BSD_license.txt   BSD License for fmlvq/fminlbfgs.m
% contrib/gmlvq/GMLVQ_classify.m   
% contrib/gmlvq/fminlbfgs.m   
% contrib/gmlvq/gmlvq.m   
% contrib/gmlvq/gmlvq_core.m   
% contrib/gmlvq/grlvq.m   
% contrib/gmlvq/grlvq_core.m   
% contrib/gmlvq/license-gpl2.txt   GNU GENERAL PUBLIC LICENSE Version 2
% 
% Contributed package RSOM
% 
% contrib/rsom/exampleDissimilarity.mat   Demo data set
% contrib/rsom/rsom_batchtrain.m   
%     contrib/rsom/rsom_demo.m   
%  contrib/rsom/rsom_lininit.m   
% contrib/rsom/rsom_randinit.m   
%     contrib/rsom/rsom_show.m   
%     contrib/rsom/rsom_umat.m   
% 
% Contributed package rneural
% 
% contrib/rneural/rneural_gas.m   quantizes the data space using the Relational Batch Neural Gas (RNG)
% 
% Other
% 
%       som/somtoolbox   this file
%     som/Contents.txt   identical to this file
%    data/winedata.mat   Wine data set (used in demos)
%       data/iris.data   IRIS data set (used in demos)
%            CHANGELOG   Changes for v2.1
%              COPYING   GNU GENERAL PUBLIC LICENSE Version 2 
%               README   Copyright notice and licence information