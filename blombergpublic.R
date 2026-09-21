# ============================================================================
# Calculate Blomberg's K for Internal Nodes
# ============================================================================
# This is a modification of the Blomberg's K to calculate phylogenetic signal 
# across the number of putative de novo events found in the inner nodes of the
# various phylogenetics trees.
# ============================================================================

library(ape)


setwd("D:/model_project/brassicaceae")
# Read NEXUS file from disk
tree <- read.tree("SpeciesTree_rooted_node_labels.txt")
class(tree)  # Should show "phylo"


setwd("D:/model_project/poaceae")
# Read NEXUS file from disk
tree <- read.tree("SpeciesTree_rooted_node_labels.txt")
class(tree) 



# Read NEXUS file
tree <- read.tree("salicaceae_tree.txt")
class(tree) 



int_nodes <- tree$node.label
de_novo_counts <- setNames(rpois(length(int_nodes), lambda=5), int_nodes)



calculate_blomberg_K_internal<- function(tree, internal_trait_values) {
  
  # Step 1: Identify internal nodes (excluding root)
  n_tips <- Ntip(tree)
  root_node <- n_tips + 1
  all_internal <- (n_tips + 1):(n_tips + Nnode(tree))
  internal_nodes <- all_internal[all_internal != root_node]
  
  
  
  # Validate input
  valid_nodes <- intersect(as.numeric(names(internal_trait_values)), internal_nodes)
  
  if(length(valid_nodes) < 2) {
    stop("Need at least 2 non-root internal nodes with trait data")
  }
  
  valid_names <- as.character(valid_nodes)
  trait_ordered <- internal_trait_values[valid_names]
  n_int <- length(valid_nodes)
  
  cat("Using", n_int, "internal nodes (excluding root)\n")
  
  # Step 2: Get distances between ALL nodes using dist.nodes()
  # This returns (Ntip + Nnode) × (Ntip + Nnode) matrix
  all_distances <- dist.nodes(tree)
  cat("Distance matrix dimensions:", dim(all_distances), "\n")
  
  # Step 3: Extract internal node submatrix
  internal_dist <- all_distances[valid_nodes, valid_nodes]
  
  # Step 4: Get heights of internal nodes
  heights <- node.depth.edgelength(tree)
  valid_heights <- heights[valid_nodes]
  
  # Step 5: Convert distances to covariances
  V_internal <- matrix(0, n_int, n_int)
  for(i in 1:n_int) {
    for(j in 1:n_int) {
      # Covariance = (height_i + height_j - distance_ij) / 2
      V_internal[i,j] <- (valid_heights[i] + valid_heights[j] - 
                            internal_dist[i,j]) / 2
    }
  }
  rownames(V_internal) <- colnames(V_internal) <- valid_names
  
  # Step 6: Ensure positive definiteness
  # Check for zero variance (root-like nodes)
  min_diag <- min(diag(V_internal))
  if(min_diag <= 0) {
    epsilon <- max(abs(V_internal)) * 1e-8 + 1e-10
    diag(V_internal) <- diag(V_internal) + epsilon
  }
  
  # Step 7: Calculate Blomberg's K using GLS
  ones <- rep(1, n_int)
  
  # Use tryCatch for matrix inversion
  inv_V <- tryCatch({
    solve(V_internal)
  }, error = function(e) {
    # If singular, add regularization
    cat("Matrix singular, adding regularization...\n")
    lambda <- 0.01
    V_reg <- V_internal + lambda * diag(diag(V_internal))
    solve(V_reg)
  })
  
  # GLS mean
  mu_hat <- as.numeric((t(ones) %*% inv_V %*% trait_ordered) / 
                         (t(ones) %*% inv_V %*% ones))
  
  # Phylogenetic MSE
  residuals_phy <- trait_ordered - mu_hat
  MSE_phy <- as.numeric((t(residuals_phy) %*% inv_V %*% residuals_phy) / 
                          (n_int - 1))
  
  # OLS MSE
  mu_ols <- mean(trait_ordered)
  residuals_ols <- trait_ordered - mu_ols
  MSE_ols <- sum(residuals_ols^2) / (n_int - 1)
  
  # Observed ratio
  observed_ratio <- MSE_ols / MSE_phy
  
  # Expected ratio under Brownian motion
  trace_V <- sum(diag(V_internal))
  mean_V <- sum(V_internal) / n_int^2
  expected_ratio <- (trace_V / n_int - mean_V) / ((n_int - 1) / n_int)
  
  K_observed <- observed_ratio / expected_ratio
  
  # Step 8: Permutation test
  n_perm <- 1000
  K_perm <- numeric(n_perm)
  
  for(i in 1:n_perm) {
    shuffled_trait <- sample(trait_ordered)
    
    mu_perm <- as.numeric((t(ones) %*% inv_V %*% shuffled_trait) / 
                            (t(ones) %*% inv_V %*% ones))
    resid_perm <- shuffled_trait - mu_perm
    MSE_perm <- as.numeric((t(resid_perm) %*% inv_V %*% resid_perm) / 
                             (n_int - 1))
    
    ratio_perm <- MSE_ols / MSE_perm
    K_perm[i] <- ratio_perm / expected_ratio
  }
  
  p_value <- (sum(K_perm >= K_observed) + 1) / (n_perm + 1)
  
  return(list(
    K = K_observed,
    p_value = p_value,
    observed_ratio = observed_ratio,
    expected_ratio = expected_ratio,
    n_internal = n_int,
    perm_dist = K_perm
  ))
}



#this calculates the internal nodes of the tree so it can be run without problems by 
#the calculate_blomberg_K_internal_corrected

n_tips <- Ntip(tree)
n_internal <- Nnode(tree)
root_node <- n_tips + 1
internal_nodes <- (n_tips + 2):(n_tips + n_internal)
names(de_novo_counts) <- as.character(internal_nodes)

#the function is applied over the tree and the de novo events per node

result <- calculate_blomberg_K_internal(tree, de_novo_counts)
result$K
result$p_value

