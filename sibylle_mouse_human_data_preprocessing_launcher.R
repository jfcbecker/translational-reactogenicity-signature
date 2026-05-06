#######################################################################
### ALL PREPROCESSING STEPS FOR THE BIOVACSAFE MOUSE AND HUMAN DATA ###
#######################################################################




### Data preprocessing - mouse data
source(file.path(this.path::here(),
                 "sibylle_mouse_data_preprocessing.R"))


### Data preprocessing - human data
source(file.path(this.path::here(),
                 "sibylle_human_data_preprocessing.R"))


### Mouse-human orthology - preamble
source(file.path(this.path::here(),
                 "sibylle_mouse_human_orthology_preamble.R"))


### Mouse-human orthology - make orthology
source(file.path(this.path::here(),
                 "sibylle_mouse_human_orthology.R"))

