.onAttach <- function(lib, pkg)
{
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  if(interactive())
  {
    #Rfiglet::figlet(message = "raedda", font = "roman")
    packageStartupMessage(
      "
                            _
                           (_ )
 _ __  _   _  _ _      ___  | |    _ _   ___   ___
( '__)( ) ( )( '_`\\  /'___) | |  /'_` )/',__)/',__)        Robust Updating
| |   | (_) || (_) )( (___  | | ( (_| |\\__, \\__, \\         Classification Rules
(_)   `\\___/'| ,__/'`\\____)(___)`\\__,_)(____/(____/
             | |
             (_)

      ", "version ", version, "\n" )
  }
  else
  { packageStartupMessage("Package 'rupclass' version ", version) }

  packageStartupMessage("Type 'citation(\"rupclass\")' for citing this R package in publications.")
  invisible()
}
