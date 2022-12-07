#Write function that calculates tip rates

GetTipRates <- function(tree = NULL,
                        Q = NULL,
                        tip.states = NULL,
                        hyper = FALSE,
                        p.mat = NULL){
  
  ### --- define inputs --- ###
  # tree: a rooted tree of class phylo
  # Q: a transition matrix
  # tip.states: a named vector of size Ntips specifying the state of each tip in terms
        # of an integer from 1 to Nstates where the names are species names that
        # match the tip labels in the provided tree. 
  # hyper: TRUE if you want the model to include a binary hyper state. FALSE if
        # you do not want the model to include a binary hyper state.
  # p.mat: a probability matrix where each column represent a discrete state and 
        #each row is a species. Values in the matrix describe the probability of
        #being in any of these states. 
 
  ### --- set default parameters --- ###
  

  ### --- perform checks --- ###
  
  #if the tree is not ultrametric, stop function and ask user to resolve
  if(is.ultrametric(tree) == FALSE){
    stop("\n Tree is not ultrametric. Please use a formal method to ultrametricize your tree. 
        If your tree is failing is.ultrametric due to rounding, you can use force.ultrametric, 
        but this does not serve as a substitute for formal rate-smoothing methods.")
  }
  
  #if tip states are null, stop function and ask user to resolve
  if(is.null(tip.states) && is.null(p.mat)){
    stop("\n Tip states are not present. Please provide tip states as integers or provide
         a probability matrix to create tip states for given data.")
  }
  
  #if tip states are not integers, stop function and ask user to resolve
  if(is.integer(tip.states)){
    stop("\n Tip states are not integers. Please provide tip states as integers.")
  }
  
  #if tip states are not present, but a probability matrix is provided, fill in 
  # tip states using the given probability matrix
  if(is.null(tip.states) && !is.null(p.mat)){
    #create vector to store tip states in
    tip.states <- c()
    #create tip states from probability matrix
    for(i in 1:nrow(p.mat)){
      tip.states[i] <- which(p.mat[i,] == 1)}
    names(tip.states) <- row.names(p.mat)
  }
  
  #loop through to store whether tip states and tree tip labels are matching for
  #each tip on the tree
 
  #empty vector to store the correct order in
  hit <- c()
  #loop through to store the correct order of the data to match the tree tips
  for(j in 1:length(tree$tip.label)){
    hit[j] <- which(tree$tip.label[j] == names(tip.states))
  }
  #reorder the data into the data frame using the correct order
  tip.states <- tip.states[hit]
  
  #check to see if the tip states and tree tip labels are in the same order
  
  print("Checking the order of tree tip labels and provided data.")
  
  order <- c()
  for(i in 1:length(tree$tip.label)){
    order[i] <- tree$tip.label[i] == names(tip.states)[i]
  }
  
  if(sum(order) == length(tree$tip.label)){
    print("Tree tip labels and provided data are in the correct order.")
  } else{
    print("Tree tip labels and provided data are not in the correct order.")
  }
  

  ### --- perform ancestral state reconstruction --- ###
  
  #reconstruct the ancestral states
  recon <- asr_mk_model(tree = tree,
                       tip_states = tip.states,
                       transition_matrix = Q,
                       Nstates = ncol(Q),
                       include_ancestral_likelihoods = T)
  
  ### --- back transform the states from 1:Nstates into their actual state --- ###
  
  # change back the state names to the actual chromosome numbers for tip rate 
  # calculation
  colnames(recon$ancestral_likelihoods) <- colnames(Q)
  
  # if there is a hyper state associated with the model, run additional code to 
  # calculate tip rates without the hyperstate labels
  if(hyper == TRUE){
    for(i in 1:ncol(recon$ancestral_likelihoods)){
      colnames(recon$ancestral_likelihoods) <- as.numeric(gsub("h", 
                                                               "", 
                                                               colnames(recon$ancestral_likelihoods)))
    }
  }
  
  #store the state names from the ASR in a vector
  states <- colnames(recon$ancestral_likelihoods)
  
  #vector to store the transformed tip states
  tip.state <- c()
  
  #loop through the integer tip states from 1:Nstates and change them back to 
  #the actual tip states supplied
  for(i in 1:length(tip.states)){
    #grabs the tip state that matches the modified tip state for the ASR
    tip.state[i] <- states[tip.states[i]]
    #ensures that the tip states are numeric 
    tip.state <- as.numeric(tip.state)
  }
  
  
  ### --- gather data for parent and tip states for tip rate calculation--- ###
  
  #create a vector to store likelihoods at each node
  est <- c()
  #loop through to store the maximum likelihood at each node
  for(i in 1:nrow(recon$ancestral_likelihoods)){
    #grabs the maximum likelihood at each node within the supplied tree
    est[i] <- which.max(as.vector(recon$ancestral_likelihoods[i,]))
    names(est) <- seq(from = 1+length(tree$tip.label), to = length(tree$tip.label)+length(est))
  }
  
  #get the node numbers for each of the tips
  tips <- tree$tip.label
  #store the node number that matches each tip label in a named vector
  nodes <- c()
  #loop through to store the tip label node numbers
  for(i in 1:length(tree$tip.label)){
    nodes[i] <- which(tree$tip.label[i] == tips)
  }
  #name the nodes
  names(nodes) <- names(tip.states)
    
  #get the edge lengths for each node
  edge.lengths <- c()
  #loop through to store the edge lengths that correspond to each node
  edge.lengths <- setNames(tree$edge.length[sapply(nodes,
                                                    function(x,y) which(y==x),y=tree$edge[,2])],names(nodes))
  
  #store the parent nodes for each tip
  nodepulls <- c()
  #loop through to store the parent node of each tip
  for(i in 1:length(tree$tip.label)){
    nodepulls[i] <- getParent(tree, nodes[i])
  }
  
  #store the state at each of the nodes
  node.state <- c()
  node.states <- c()
  #loop through to store the state corresponding to each node
  for(i in 1:length(nodepulls)){
    node.state[i] <- which(names(est) == nodepulls[i])
    #pull the 
    node.states[i] <- est[node.state[i]]
    #pull the actual chromosome number value for each node
    node.states[i] <- states[node.states[i]]
    #ensures that the node states are numeric 
    node.states <- as.numeric(node.states)
  }
  
  #store the tip rate for each of the tips
  tip.changes <- c()
  #loop through to calculate the tip rate for each tip on the tree, given 
  #the node and tip state information
  for(i in 1:length(tip.state)){
    #calculate the tip rates
    tip.changes[i] <- abs((tip.state[i]-node.states[i]))/edge.lengths[i]
  }
  names(tip.changes) <- tree$tip.label
return(tip.changes)
}










