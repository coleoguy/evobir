ViewEvo <- function(simulation){
  if(simulation == 'wf.model') shiny::runApp(system.file("wf.model", package='evobiR'))
  if(simulation == 'bd.model') shiny::runApp(system.file("bd.model", package='evobiR'))
  if(simulation == 'dist.model') shiny::runApp(system.file("dist.model", package='evobiR'))
  if(simulation == 'bm.model') shiny::runApp(system.file("bm.model", package='evobiR'))
  if(simulation == 'clumping.model') shiny::runApp(system.file("clumping.model", package='evobiR'))
}