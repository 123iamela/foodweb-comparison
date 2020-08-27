# Select the Weddell Sea Food Web 
#
require(dplyr)

tbl <- read.delim("Data/bodysizes_2008.txt",stringsAsFactors = FALSE)

tbl1 <- read.csv("Data/SothernOceanDiets.csv",stringsAsFactors = FALSE)

tbl2 <-filter(tbl1,  grepl("Amphipoda",PREY_NAME))
tbl2 <-filter(tbl1,  grepl("Cetacea",PREDATOR_COMMON_NAME))
tbl2 <-filter(tbl1,  grepl("Weddell",LOCATION))

names(tbl1)
unique(tbl1$LOCATION)

names(tbl)

tbl <-filter(tbl,  grepl("Weddel",Geographic.location))

unique(tbl$Geographic.location)

tbl <- select(tbl, Link.ID, Taxonomy.consumer,Lifestage.consumer,Taxonomy.resource,Lifestage...resource)

unique(tbl$Lifestage.consumer)
unique(tbl$Lifestage...resource)
unique(tbl$Metabolic.category.resource)
tbl <- filter(tbl, Lifestage...resource!="pups")

tbl <- select(tbl, Link.ID, Taxonomy.consumer,Taxonomy.resource)


write.csv(tbl,"Weddel.csv",col.names = TRUE)


#
# Select marine 
# 

tbl <- read.delim("Data/bodysizes_2008.txt",stringsAsFactors = FALSE)

names(tbl)


unique(tbl$Geographic.location)

mar <-filter(tbl,  grepl("marine",General.habitat))

unique(mar$Geographic.location)
