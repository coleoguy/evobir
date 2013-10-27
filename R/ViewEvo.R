ViewEvo <- function(x){
  if(x == 'wf.model') shiny::runApp(system.file("wf.model", package='evobiR'))
  if(x == 'bd.model') shiny::runApp(system.file("bd.model", package='evobiR'))
}