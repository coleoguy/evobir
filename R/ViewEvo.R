ViewEvo <- function(simulation){
  if(simulation == 'wf.model') runApp(system.file("wf.model", package='evobiR'))
  if(simulation == 'bd.model') runApp(system.file("bd.model", package='evobiR'))
  if(simulation == 'dist.model') runApp(system.file("dist.model", package='evobiR'))
  if(simulation == 'bm.model') runApp(system.file("bm.model", package='evobiR'))
  if(simulation == 'clumping.model') runApp(system.file("clumping.model", package='evobiR'))
  if(simulation == 'bm.tree.model') runApp(system.file("bm.tree.model", package='evobiR'))
}