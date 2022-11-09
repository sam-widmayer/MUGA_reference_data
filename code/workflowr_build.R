library(workflowr)

# Initializing workflowr infrastructure
wflow_start(directory = "~/Documents/MUGA_reference_data/", 
            existing = TRUE) # adds workflowr files to existing directory

# Preview your changes
wflow_build()
