# Test mangal

require(rmangal)
require(plyr)
require(dplyr)
netdb <- mangalapi()
netdb
whatIs(netdb, 'taxa')
whatIs(netdb, 'network')
whatIs(netdb, 'dataset')
data_set <-   
nets <-ldply(listNetwork(netdb), summarize, id = id, description = description, name = name, n_int = length(interactions), 
             papers = )

