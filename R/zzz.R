.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\nThank you for using evobiR!\n")
  packageStartupMessage("\nTo acknowledge our work, please cite the package: \n")
  packageStartupMessage(" Michelle M. Jonika, Maximos Chin, Nathan Anderson, Richard H. Adams, Jeffery P. Demuth, and Heath Blackmon. (2025).")
  packageStartupMessage(" EvobiR: Tools for comparative analyses and teaching evolutionary biology.")
  packageStartupMessage(" coleoguy/evobiR: EvobiR version 2.2 \n")
}
