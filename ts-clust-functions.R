

getBestKPerCVI <- function(cvis){
  
  #to_be_minimized <- c("DB", "DBstar", "COP")
  to_be_minimized <- c("DB", "DBstar", "COP", "VI", "K", "T")
  
  #to_be_maximized <- c("Sil", "SF", "CH", "D")
  to_be_maximized <- c("Sil", "SF", "CH", "D", "MPC", "SC", "PBMF")
  
  k_range <- as.numeric( gsub("k_", "", colnames(cvis), fixed = TRUE) )
  
  
  best_idx_per_cvi <- c(
    # apply( cvis[to_be_minimized,], MARGIN = 1, function(row) which.min(row[row > 0]) ),
    # apply( cvis[to_be_maximized,], MARGIN = 1, function(row) which.max(row) )
    apply( cvis[rownames(cvis) %in% to_be_minimized,], MARGIN = 1, 
           function(row) which.min(row[row > 0]) ),
    apply( cvis[rownames(cvis) %in% to_be_maximized,], MARGIN = 1, 
           function(row) which.max(row) )
  )
  
  best_k_per_cvi <- sapply( best_idx_per_cvi, function(i) k_range[i] )
  
  votes_for_k_overall <- table(best_k_per_cvi)
  best_k_overall <- as.numeric( names( which( votes_for_k_overall == max(votes_for_k_overall)) ) )
  
  if ( length(best_k_overall) > 1 ) {
    
    best_idx_per_best_ks <- c(
      apply( cvis[rownames(cvis) %in% to_be_minimized,paste0("k_", best_k_overall)], MARGIN = 1, 
             function(row) which.min(row) ),
      apply( cvis[rownames(cvis) %in% to_be_maximized,paste0("k_", best_k_overall)], MARGIN = 1, 
             function(row) which.max(row) )  
    )
    
    best_selected_k_per_cvi <- sapply( best_idx_per_best_ks, function(i) best_k_overall[i] )
    
    votes_for_best_ks <- table(best_selected_k_per_cvi)
    best_k <- as.numeric( names( which( votes_for_best_ks == max(votes_for_best_ks)) ) )
    
  } else {
    
    best_selected_k_per_cvi <- NULL
    best_k <- best_k_overall
    
  }
  
  
  list(
    best_k_per_cvi = best_k_per_cvi,
    best_selected_k_per_cvi = best_selected_k_per_cvi,
    best_k = best_k)
}


tsClusterings <- function(data, k_range, trace = FALSE, n_rep = 5, seed) {
  
  elaboration_start_date <- now()
  
  
  pc_dtw_dba <- tsclust(data, k = k_range,
                        distance = "dtw_basic", centroid = "dba",
                        trace = trace, seed = seed,
                        norm = "L2",
                        #control = partitional_control(nrep = n_rep),
                        args = tsclust_args(cent = list(trace = trace)) )
  
  pc_dtw_pam <- tsclust(data, k = k_range,
                        distance = "dtw_basic", centroid = "pam",
                        trace = trace, seed = seed,
                        norm = "L2",
                        #control = partitional_control(nrep = n_rep),
                        args = tsclust_args(cent = list(trace = trace)) )
  
  
  pc_gak_dba <- tsclust(data, k = k_range,
                        distance = "gak", centroid = "dba",
                        trace = trace, seed = seed,
                        norm = "L2",
                        #control = partitional_control(nrep = n_rep),
                        args = tsclust_args(cent = list(trace = trace)) )
  
  
  pc_gak_pam <- tsclust(data, k = k_range,
                        distance = "gak", centroid = "pam",
                        trace = trace, seed = seed,
                        norm = "L2",
                        #control = partitional_control(nrep = n_rep),
                        args = tsclust_args(cent = list(trace = trace)) )
  
  
  pc_softdtw_pam <- tsclust(data, k = k_range,
                            distance = "soft-dtw", centroid = "pam",
                            trace = trace, seed = seed,
                            norm = "L2",
                            #control = partitional_control(nrep = n_rep),
                            args = tsclust_args(cent = list(trace = trace)) )
  
  pc_softdtw_sdtw <- tsclust(data, k = k_range,
                             distance = "soft-dtw", centroid = "sdtw_cent",
                             trace = trace, seed = seed,
                             norm = "L2",
                             #control = partitional_control(nrep = n_rep),
                             args = tsclust_args(cent = list(trace = trace)) )
  
  
  # names(pc_dtw_dba_rep) <- paste0("rep_", 1L:n_rep)
  # names(pc_dtw_pam_rep) <- paste0("rep_", 1L:n_rep)
  # names(pc_gak_dba_rep) <- paste0("rep_", 1L:n_rep)
  # names(pc_gak_pam_rep) <- paste0("rep_", 1L:n_rep)
  # names(pc_softdtw_pam_rep) <- paste0("rep_", 1L:n_rep)
  
  # pc_dtw_dba_ens <- cl_ensemble(list = pc_dtw_dba_rep)
  # pc_dtw_pam_ens <- cl_ensemble(list = pc_dtw_pam_rep)
  # pc_gak_dba_ens <- cl_ensemble(list = pc_gak_dba_rep)
  # pc_gak_pam_ens <- cl_ensemble(list = pc_gak_pam_rep)
  # pc_softdtw_pam_ens <- cl_ensemble(list = pc_softdtw_pam_rep)
  
  names(pc_dtw_dba) <- paste0("k_", k_range)
  names(pc_dtw_pam) <- paste0("k_", k_range)
  names(pc_gak_dba) <- paste0("k_", k_range)
  names(pc_gak_pam) <- paste0("k_", k_range)
  names(pc_softdtw_pam) <- paste0("k_", k_range)
  names(pc_softdtw_sdtw) <- paste0("k_", k_range)
  
  
  ##  Evaluations
  pc_dtw_dba_cvis <- sapply(pc_dtw_dba, cvi)
  pc_dtw_pam_cvis <- sapply(pc_dtw_pam, cvi)
  pc_gak_dba_cvis <- sapply(pc_gak_dba, cvi)
  pc_gak_pam_cvis <- sapply(pc_gak_pam, cvi)
  pc_softdtw_pam_cvis <- sapply(pc_softdtw_pam, cvi)
  pc_softdtw_sdtw_cvis <- sapply(pc_softdtw_sdtw, cvi)
  
  
  pc_dtw_dba_best_k <- getBestKPerCVI( pc_dtw_dba_cvis )
  pc_dtw_pam_best_k <- getBestKPerCVI( pc_dtw_pam_cvis )
  pc_gak_dba_best_k <- getBestKPerCVI( pc_gak_dba_cvis )
  pc_gak_pam_best_k <- getBestKPerCVI( pc_gak_pam_cvis )
  pc_softdtw_pam_best_k <- getBestKPerCVI( pc_softdtw_pam_cvis )
  pc_softdtw_sdtw_best_k <- getBestKPerCVI( pc_softdtw_sdtw_cvis )
  
  
  compare_best_k <- cbind(
    pc_dtw_dba = pc_dtw_dba_best_k$best_k_per_cvi,
    pc_dtw_pam = as.integer(pc_dtw_pam_best_k$best_k_per_cvi),
    pc_gak_dba = pc_gak_dba_best_k$best_k_per_cvi,
    pc_gak_pam = as.integer(pc_gak_pam_best_k$best_k_per_cvi),
    pc_softdtw_pam = as.integer(pc_softdtw_pam_best_k$best_k_per_cvi),
    pc_softdtw_sdtw = as.integer(pc_softdtw_sdtw_best_k$best_k_per_cvi)
  )
  
  elaboration_end_date <- now()
  
  elaboration_duration <- as.duration( elaboration_start_date %--% elaboration_end_date )
  
  
  # Output
  
  list(
    clusterings = list( 
      dtw_dba = list( tsclustering = pc_dtw_dba,
                      cvis = pc_dtw_dba_cvis,
                      best_k = pc_dtw_dba_best_k ),
      
      dtw_pam = list( tsclustering = pc_dtw_pam,
                      cvis = pc_dtw_pam_cvis,
                      best_k = pc_dtw_pam_best_k ),
      
      gak_dba = list( tsclustering = pc_gak_dba,
                      cvis = pc_gak_dba_cvis,
                      best_k = pc_gak_dba_best_k ),
      
      gak_pam = list( tsclustering = pc_gak_pam,
                      cvis = pc_gak_pam_cvis,
                      best_k = pc_gak_pam ),
      
      softdtw_pam = list( tsclustering = pc_softdtw_pam,
                          cvis = pc_softdtw_pam_cvis,
                          best_k = pc_softdtw_pam_best_k ),
      
      softdtw_sdtw = list( tsclustering = pc_softdtw_sdtw,
                           cvis = pc_softdtw_sdtw_cvis,
                           best_k = pc_softdtw_sdtw_best_k ) 
    ),
    
    compare_best_k = compare_best_k,
    compare_best_k_df = as_tibble(compare_best_k) %>%
      bind_cols( CVI = rownames(compare_best_k) ),
    
    most_voted_k = sort( table(compare_best_k), decreasing = TRUE ),
    most_voted_k_df = as_tibble( sort( table(compare_best_k), decreasing = TRUE ) ),
    
    duration = elaboration_duration
  )
}


bestKGreaterThanKMin <- function(clusterings, k_min) {
  ks_greater_than_k_min <- clusterings$most_voted_k[ as.numeric(names(clusterings$most_voted_k)) >= k_min ]
  
  best_ks_greater_than_k_min <- as.numeric(
    names( which(ks_greater_than_k_min == max(ks_greater_than_k_min)) )
  )
  
  names(best_ks_greater_than_k_min) <- as.character(best_ks_greater_than_k_min)
  
  most_voted_k_greather_then_k_min <-
    map_df( best_ks_greater_than_k_min,
            function(k) sum( map_dbl(clusterings$clusterings,
                                     ~ sum(.$best_k$best_selected_k_per_cvi == k)) ) )
  
  selected_k <- as.numeric(
    names( sort(as_vector(most_voted_k_greather_then_k_min), decreasing = TRUE)[1] )
  )
  
  selected_k
}


extractCentroids <- function(data, clusterings, selected_k, n_rep, seed = NULL){
  
  global_weighted_av_dist <- function(clus_info_df) {
    clus_info_df %>% 
      mutate( av_dist_tot = sum(av_dist) ) %>% 
      mutate( av_dist_perc = av_dist / av_dist_tot ) %>% 
      mutate( av_dist_weighted = size * av_dist_perc ) %>% 
      mutate( global_weighted_av_dist = sum(av_dist_weighted) / sum(size) ) %>% 
      distinct( global_weighted_av_dist ) %>% 
      pull()  
  }
  
  # What's the method that gives the selected k for the most number of times
  # in its CVIs?
  selected_algorithm_counts <- sort( colSums( clusterings$compare_best_k == selected_k ), decreasing = TRUE )
  
  selected_algorithms <- names( which( selected_algorithm_counts == max(selected_algorithm_counts) ) )
  
  
  global_weighted_av_dist_per_alg <- c(
    
    global_weighted_av_dist(
      clusterings$clusterings$dtw_dba$tsclustering[[paste0("k_", selected_k)]]@clusinfo ),
    
    global_weighted_av_dist( 
      clusterings$clusterings$dtw_pam$tsclustering[[paste0("k_", selected_k)]]@clusinfo ),
    
    global_weighted_av_dist( 
      clusterings$clusterings$gak_dba$tsclustering[[paste0("k_", selected_k)]]@clusinfo ),
    
    global_weighted_av_dist( 
      clusterings$clusterings$gak_pam$tsclustering[[paste0("k_", selected_k)]]@clusinfo ),
    
    global_weighted_av_dist( 
      clusterings$clusterings$softdtw_pam$tsclustering[[paste0("k_", selected_k)]]@clusinfo ),
    
    global_weighted_av_dist( 
      clusterings$clusterings$softdtw_sdtw$tsclustering[[paste0("k_", selected_k)]]@clusinfo )
    
  )
  
  names(global_weighted_av_dist_per_alg) <- c("pc_dtw_dba", "pc_dtw_pam", "pc_gak_dba",
                                              "pc_gak_pam", "pc_softdtw_pam", "pc_softdtw_sdtw")
  
  
  selected_global_dist <- global_weighted_av_dist_per_alg[selected_algorithms]
  
  pc_ens_centroid_list <- list()
  
  for (alg in names(selected_global_dist)) {
    pc_rep_str <- case_when(
      
      alg == "pc_dtw_pam" ~ paste0("tsclust(data, k = selected_k,
                                             distance = \"dtw_basic\", centroid = \"pam\",
                                             trace = FALSE, seed = ", seed, ",
                                             norm = \"L2\",
                                             control = partitional_control(nrep = n_rep) )"),
      
      alg == "pc_dtw_dba" ~ paste0("tsclust(data, k = selected_k,
                                             distance = \"dtw_basic\", centroid = \"dba\",
                                             trace = FALSE, seed = ", seed, ",
                                             norm = \"L2\",
                                             control = partitional_control(nrep = n_rep) )"),
      
      alg == "pc_gak_pam" ~ paste0("tsclust(data, k = selected_k,
                                             distance = \"gak\", centroid = \"pam\",
                                             trace = FALSE, seed = ", seed, ",
                                             norm = \"L2\",
                                             control = partitional_control(nrep = n_rep) )"),
      
      alg == "pc_gak_dba" ~ paste0("tsclust(data, k = selected_k,
                                             distance = \"gak\", centroid = \"dba\",
                                             trace = FALSE, seed = ", seed, ",
                                             norm = \"L2\",
                                             control = partitional_control(nrep = n_rep) )"),
      
      alg == "pc_softdtw_pam" ~ paste0("tsclust(data, k = selected_k,
                                                 distance = \"soft-dtw\", centroid = \"pam\",
                                                 trace = FALSE, seed = ", seed, ",
                                                 norm = \"L2\",
                                                 control = partitional_control( nrep = n_rep) )"),
      
      alg == "pc_softdtw_sdtw" ~ paste0("tsclust(data, k = selected_k,
                                                 distance = \"soft-dtw\", centroid = \"sdtw_cent\",
                                                 trace = FALSE, seed = ", seed, ",
                                                 norm = \"L2\",
                                                 control = partitional_control( nrep = n_rep) )")
    )
    
    pc_rep <- eval( parse(text = pc_rep_str) )
    
    
    names(pc_rep) <- paste0("rep_", 1L:n_rep)
    
    pc_ens <- cl_ensemble(list = pc_rep)
    
    pc_ens_centroid <- cl_medoid(pc_ens)
    
    pc_ens_centroid_list[alg] <- pc_ens_centroid
  }
  
  
  list(
    selected_global_distances = selected_global_dist,
    pc_ens_list = pc_ens_centroid_list,
    selected_pc_ens = pc_ens_centroid_list[[ names(which.min(selected_global_dist)) ]]
  )
}

extractClusteringResultsCSV <- function(medoids, file_path = "data_clusterized.csv", data, extract_CSV = TRUE) {
  
  selected_k <- medoids@k
  
  cluster_column_name <- paste0("cluster_id_k_", selected_k)
  medoid_column_name <- paste0("is_medoid_k_", selected_k)
  
  clustering_result <- tibble(
    timeseries_key_column_name = names(medoids@datalist),
    cluster_id = medoids@cluster,
    is_medoid = ifelse( as.vector(medoids@cldist) == 0, 1, 0 )
  )
  
  timeseries_sorting <- tibble(
    timeseries_key_column_name = sort( names(medoids@datalist) ),
    interactive_sorting = length(medoids@datalist)
  )
  
  clustering_result_sorted <- clustering_result %>%
    inner_join(timeseries_sorting, by = "timeseries_key_column_name")
  
  
  data_clusterized_tbl <- data %>%
    inner_join(clustering_result, by = setNames(nm=names(data)[1], "timeseries_key_column_name")) %>%
    rename( !!cluster_column_name := cluster_id,
            !!medoid_column_name  := is_medoid )
  
  # Export the clustering result in csv
  if (extract_CSV){
    data_clusterized_tbl %>%
      write_csv(path = file_path)
  }
  
  data_clusterized_tbl
}

