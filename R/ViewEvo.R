ViewEvo <- function(simulation){
  mods <- c('wf.model', 'bd.model', 
            'dist.model','bm.model', 
            'clumping.model', 'bm.tree.model',
            'treeviz')
  if(simulation == 'wf.model') runApp(system.file("wf.model", package='evobiR'))
  if(simulation == 'bd.model') runApp(system.file("bd.model", package='evobiR'))
  if(simulation == 'dist.model') runApp(system.file("dist.model", package='evobiR'))
  if(simulation == 'bm.model') runApp(system.file("bm.model", package='evobiR'))
  if(simulation == 'clumping.model') runApp(system.file("clumping.model", package='evobiR'))
  if(simulation == 'bm.tree.model') runApp(system.file("bm.tree.model", package='evobiR'))
  if(simulation == 'treeviz') runApp(system.file("treeviz", package='evobiR'))
  if(!simulation %in% mods){
    cat("The specified app was not found in evobiR.  Available apps are:")
    cat(mods)
  }
}


