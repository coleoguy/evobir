#### FIX SIMMAP  ####

#Arguments: hists(simmap or multisimmap object to be fixed)
#           tips(two column dataframe with first column listing species name as
#                shown in tree tips and second column listing simulation state
#                for each species)
#           transition.matrix(symmetrical matrix describing possible transitions
#                             between states)
#
#Output: simmap or multisimmap object with failing edges fixed (time spent in 
#        each transitional state is assumed to be equal), maps and mapped.edges
#        adjusted accordingly

fix.simmap <- function(hists, tips, transition.matrix){
  
  #### REARRANGE IF SINGLE MAP ####
  if("simmap"  %in% class(hists)){
    hist <- hists
    hists <- list()
    hists[[1]] <- hist
  }
  #### CATALOG FAILED EDGES ####
  
  # Create list of 100 lists to store failures
  fail.list <- rep(list(c()), length(hists))
  
  #Loop through 100 simulations and store failed tips
  for (i in 1:length(hists)){
    for (j in 1:nrow(hists[[i]]$edge)){
      if("fail" %in% names(hists[[i]]$maps[[j]])){
        fail.list[[i]] <- c(fail.list[[i]] ,j)
      } else {
        names(hists[[i]]$maps[[j]])
      }
    }
  }
  
  #### LOOP TO FIX ####
  
  for (i in 1:length(hists)){
    
    #Check if any edges failed
    if(length(fail.list[[i]]) == 0){
      
      #Go to next simulation
      next
      
    }
    
    
    # catalog the simulation states at each node and tip of the maps
    tip.states <- c()
    
    for(j in 1:length(hists[[i]]$tip.label)){
      tip.leading.edge <- which(hists[[i]]$edge[,2] == j)
      tip.state <- names(hists[[i]]$maps[[tip.leading.edge]]
                         [length(hists[[i]]$maps[[tip.leading.edge]])])
      tip.states <- c(tip.states,tip.state)
    }
    
    names(tip.states) <- hists[[i]]$tip.label
    
    # cataloging which tips failed, the taxa on those tips, the branches leading to
    # those tips, and the simulation state from which the transition failed
    fail.taxa <- tips[which(tips[,1] %in%
                              names(which(tip.states == "fail"))), ]
    
    #Check if any tips failed
    if(nrow(fail.taxa) != 0 ){
      #add tip numbers to fail.taxa
      for(j in 1:nrow(fail.taxa)){
        fail.taxa[j,3] <- which(hists[[i]]$tip.label %in% 
                                  fail.taxa[j,1])
      }
      
      fail.tips.names <- which(tip.states == "fail")
      fail.tips.numbers <- matrix(hists[[i]]$edge[fail.list[[i]], ],
                                  ncol = 2)
      fail.tips.edges <- fail.list[[i]][which(fail.tips.numbers[,2]
                                              <= length(hists[[i]]$tip.label))]
      fail.tips.numbers <- matrix(fail.tips.numbers[which(fail.tips.numbers[,2]
                                                          <= length(hists[[i]]$tip.label)),],
                                  ncol=2)
      fail.tips.start <- fail.tips.numbers[, 1]
      fail.tips.leading.edge <- c()
      for(j in 1:length(fail.tips.start)){
        fail.tips.leading.edge <- c(fail.tips.leading.edge,which(hists[[i]]$edge[,2] == 
                                                                   fail.tips.start[j]))
      }
      fail.tips.leading.edge <- unique(fail.tips.leading.edge)
      fail.tips.leading.maps <- hists[[i]]$maps[fail.tips.leading.edge]
      fail.tips.start.states <- c()
      
      for(j in 1:length(fail.tips.leading.maps)){
        
        fail.tips.start.states[j] <- names(fail.tips.leading.maps[[j]])[
          length(fail.tips.leading.maps[[j]])]
      }
      
      names(fail.tips.start.states) <- unique(fail.tips.start)
      
      # cataloging the simulation state that the failed tips should have resolved to
      # ????? using only the first data entry: some taxa have data for multiple karyotypes ?????
      fail.tips.end.states <- fail.taxa[,2]
      
      names(fail.tips.end.states) <- fail.taxa[,3]
      
      # cataloging the length of the branches leading to the failed tips
      branch.lengths <- c()
      for (k in 1:nrow(fail.tips.numbers)){
        branch.lengths[k] <- hists[[i]]$edge.length[
          which(fail.tips.numbers[k, 1] == hists[[i]]$edge[ , 1]
                & fail.tips.numbers[k, 2] == hists[[i]]$edge[ , 2])]
      }
      branches <- cbind(fail.tips.numbers, branch.lengths)
    } else {
      
      #Set to blank
      fail.tips.edges <- NA
      
    }
    
    # cataloging which internal nodes failed, the branches leading to those nodes,
    # and the simulation states from which the transition failed
    fail.internal.edges <- fail.list[[i]][which(!fail.list[[i]] %in% 
                                                  fail.tips.edges)]
    
    if(length(fail.internal.edges) != 0){
      fail.internal <- as.numeric(hists[[i]]$edge[fail.internal.edges,2])
      fail.internal.start <- hists[[i]]$edge[
        which(hists[[i]]$edge[,2] %in% fail.internal), 1]
      fail.internal.leading.edge <- c()
      fail.internal.leading.maps <- list()
      for(j in 1:length(fail.internal.start)){
        if(fail.internal.start[j] != (length(hists[[i]]$tip.label) + 1)){
          fail.internal.leading.edge <- c(fail.internal.leading.edge,
                                          which(hists[[i]]$edge[,2] == 
                                                  fail.internal.start[j]))
          fail.internal.leading.maps[j] <- hists[[i]]$maps[which(hists[[i]]$edge[,2] == fail.internal.start[j])]
        } else {
          fail.internal.leading.edge <- c(fail.internal.leading.edge,
                                          0)
          sibling.leading.edge <- which(hists[[i]]$edge[,1] == hists[[i]]$edge[fail.internal.edges[j],1] &
                                          hists[[i]]$edge[,2] != hists[[i]]$edge[fail.internal.edges[j],2])
          fail.internal.leading.maps[[j]] <- hists[[i]]$maps[[sibling.leading.edge]][1]
        }
      }
      
      fail.internal.start.states <- c()
      
      for(j in 1:length(fail.internal.leading.maps)){
        
        fail.internal.start.states[j] <- names(fail.internal.leading.maps[[j]])[
          length(fail.internal.leading.maps[[j]])]
      }
      
      names(fail.internal.start.states) <- fail.internal.start
      if(TRUE %in% duplicated(names(fail.internal.start.states))){
        fail.internal.start.states <- fail.internal.start.states[-which(duplicated(names(fail.internal.start.states)))]
      }
      
      fail.internal.end.states <- rep("fail", length(fail.internal))
      names(fail.internal.end.states) <- fail.internal
      
      #Get failed.internal.numbers
      fail.internal.numbers <- matrix(hists[[i]]$edge[fail.internal.edges,],
                                      ncol=2)
      
      #Get failed.internal branch lengths and combine
      branches.internal <- cbind(fail.internal.numbers,
                                 hists[[i]]$edge.length[fail.internal.edges])
      if(exists("branches") &&
         ncol(branches) == 3){
        #bind failed branches to failed tips 
        branches <- rbind(branches,
                          branches.internal)
      } else {
        
        #reassign branches.internal to branches
        branches <- branches.internal
        
      }
    }
    
    #columns for extra length, edge id, start state, end state
    branches <- cbind(branches,c(rep(NA,nrow(branches))),
                      c(rep(NA,nrow(branches))),
                      c(rep(NA,nrow(branches))),
                      c(rep(NA,nrow(branches))),
                      c(rep(NA,nrow(branches))))
    
    colnames(branches) <- c("start",
                            "end",
                            "edge.length",
                            "extra.edge.length",
                            "edge",
                            "extra.edge",
                            "start.state",
                            "end.state")
    
    #add edge info
    if(length(fail.internal.edges) == 0){
      branches[,5] <- fail.tips.edges
    } else if(nrow(fail.taxa) == 0){
      branches[,5] <- fail.internal.edges
    } else {
      branches[,5] <- c(fail.tips.edges,fail.internal.edges)
    }
    
    # adding length of internal failed branches connected by internal node to failed tip
    # catalog known start and end states for failed branches for transition inference
    # standardizing branch length
    for (l in 1:nrow(branches)){
      
      #Determine start states
      
      #Check if leading edge also failed
      if (branches[l,1] %in% branches[,2]){
        
        #Retreive sibling edge (same start node, different end node
        sibling.edge <- which(hists[[i]]$edge[,1] %in% branches[l,1] &
                                !hists[[i]]$edge[,2] %in% branches[l,2])
        
        sibling.map <- hists[[i]]$maps[sibling.edge]
        
        #Check if sibling tip also failed
        if(names(sibling.map[[1]])[1] == "fail"){
          
          #New start from start node of leading edge
          start.1 <- names(fail.internal.start.states[
            names(fail.internal.end.states) == branches[l, 1]])
          
          #Extra length from leading edge
          extra.length <- hists[[i]]$edge.length[
            which(start.1 == hists[[i]]$edge[ , 1]
                  & branches[l, 1] == hists[[i]]$edge[ , 2])]
          
          branches[l,4] <- extra.length
          
          #document failed leading edge
          branches[l,6] <- which(start.1 == hists[[i]]$edge[ , 1]
                                 & branches[l, 1] == hists[[i]]$edge[ , 2])
          
          #Start state from leading edge
          branches[l, 7]  <- fail.internal.start.states[
            names(fail.internal.end.states) == branches[l, 1]]
        } else {
          
          #Start state from sibling edge
          branches[l,7] <- names(sibling.map[[1]])[1]
        }
      } else {
        
        #Check if tip or internal
        if(as.numeric(branches[l,2]) <= length(hists[[i]]$tip.label)){
          
          #Retrieve start from fail.tips.start.states
          branches[l, 7] <-
            fail.tips.start.states[names(fail.tips.start.states) == branches[l, 1]]
          
        } else {
          #Retrieve start from fail.internal.start.states
          branches[l, 7] <-
            fail.internal.start.states[names(fail.internal.start.states) == branches[l, 1]]
        }
      }
      
      #Determine end states
      
      #Check if tip or internal
      if(as.numeric(branches[l,2]) <= length(hists[[i]]$tip.label)){
        
        #Retreive end state from fail.tips.end.state 
        branches[l, 8] <-
          fail.tips.end.states[names(fail.tips.end.states) == branches[l, 2]]
      } else {
        
        #Check if trailing tip failed
        if(branches[l,2] %in% branches[,1]){
          
          #Test for descendent which didn't fail
          descendent.sibling.edge <- which(hists[[i]]$edge[,1] %in% branches[l,2] &
                                             !hists[[i]]$edge[,2] %in% branches[,2])
          
          #Check if non-failing descendent exists
          if(length(descendent.sibling.edge) == 0){
            
            #Set to fail (already added under "extra.edge" section)
            branches[l,8] <- NA
          } else {
            
            #Retreive map from edge
            descendent.sibling.map <- hists[[i]]$maps[descendent.sibling.edge]
            
            branches[l,8] <- names(descendent.sibling.map[[1]])[1]
          }
        } else {
          
          #Retreive one of two trailing edges (doesn't matter which one as 
          #both should start in same state)
          trailing <- which(hists[[i]]$edge[,1] %in% branches[l,2])[1]
          
          trailing.maps <- hists[[i]]$maps[trailing]
          
          branches[l,8] <- names(trailing.maps[[1]][1])
        }
      }
    }
    
    #Remove any branches with failing end states (redundant since they are already
    #included elsewhere)
    if(NA %in% branches[,8]){
      branches <- branches[-which(is.na(branches[,8])),]
    }
    
    #Change to matrix
    branches <- as.matrix(branches)
    
    #Retreive possible paths
    graph.paths <- graph_from_data_frame(
      which(transition.matrix != 0, arr.ind = TRUE))
    
    #Find paths and add to maps
    for(j in 1:nrow(branches)){
      
      # use graph/network theory to find the shortest path between simulation states
      found.path <- shortest_paths(graph.paths, branches[j, 7], branches[j, 8])[[1]][[1]]
      
      #Distinguish between connected and isolated failed tips
      if(is.na(branches[j,6])){
        
        #Calculate time spent in each state
        fail.time <- as.numeric(branches[j, 3])/length(found.path)
        
        #Add to maps
        hists[[i]]$maps[[as.numeric(branches[j,5])]] <- rep(fail.time,
                                                            length(found.path))
        
        #Name with states
        names(hists[[i]]$maps[[as.numeric(branches[j,5])]]) <- names(found.path)
        
      } else {
        
        #Calculate time spent in each state
        fail.time <- (as.numeric(branches[j, 3]) + as.numeric(branches[j,4]))/
          length(found.path)
        
        #Sum time and divide between edges
        leading.multiple <- as.integer(as.numeric(branches[j,4])/fail.time)
        trailing.multiple <- as.integer(as.numeric(branches[j,3])/fail.time)
        leading.leftover <- as.numeric(branches[j,4]) - (leading.multiple * fail.time)
        trailing.leftover <- as.numeric(branches[j,3]) - (trailing.multiple * fail.time)
        
        #Add to leading map
        hists[[i]]$maps[[as.numeric(branches[j,6])]] <- c(rep(fail.time,
                                                              leading.multiple),
                                                          leading.leftover)
        
        names(hists[[i]]$maps[[as.numeric(branches[j,6])]]) <- 
          names(found.path)[1:(leading.multiple + 1)]
        
        #Add to trailing map
        hists[[i]]$maps[[as.numeric(branches[j,5])]] <- c(trailing.leftover,
                                                          rep(fail.time,
                                                              trailing.multiple))
        
        names(hists[[i]]$maps[[as.numeric(branches[j,5])]]) <- 
          names(found.path)[(length(found.path) - 
                               trailing.multiple):length(found.path)]
        
        #Check if sibling tip failed
        if(sum(as.numeric(branches[,6]) == branches[j,6],na.rm = T) >= 1){
          
          #identify index associated with other failed sibling
          repeat.index <- which(branches[,6] == branches[j,6] &
                                  branches[,5] != branches[j,5])
          
          #Change extra.edge.length and extra.edge to NA
          branches[repeat.index,4] <- NA
          branches[repeat.index,6] <- NA
          
          #Change start state
          branches[repeat.index,7] <- names(hists[[i]]$maps[[as.numeric(
            branches[j,5])]])[1]
          
        }
      }
    }
    
    #Make changes to mapped edges
    for(j in 1:nrow(branches)){
      
      #retreive edge
      edge <- as.numeric(branches[j,5])
      
      #Loop through states
      for(k in 1:length(hists[[i]]$maps[[edge]])){
        
        #Retreive state
        state <- names(hists[[i]]$maps[[edge]])[k]
        
        if(state %in% colnames(hists[[i]]$mapped.edge)){
          #Sub into matrix
          hists[[i]]$mapped.edge[edge,state] <- hists[[i]]$maps[[edge]][k]
          
        } else {
          
          #Name
          colnames.new <- c(colnames(hists[[i]]$mapped.edge),
                            state)
          
          #Add column
          hists[[i]]$mapped.edge <- cbind(hists[[i]]$mapped.edge,c(0))
          
          #Rename
          colnames(hists[[i]]$mapped.edge) <- colnames.new
          
          #Sort to reorder
          hists[[i]]$mapped.edge <- hists[[i]]$mapped.edge[,sort(colnames(
            hists[[i]]$mapped.edge))]
          
          #Sub into matrix
          hists[[i]]$mapped.edge[edge,state] <- hists[[i]]$maps[[edge]][k]
          
        }
      }
      
      if("fail" %in% colnames(hists[[i]]$mapped.edge)){
        hists[[i]]$mapped.edge[edge,"fail"] <- 0
      }
    }
    
  }
  
  #### RETURN FIXED OBJECT ####
  if(length(hists) == 1){
    hists <- hists[[1]]
  }
  return(hists)
}

