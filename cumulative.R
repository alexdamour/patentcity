source("loadGraphs.R")

bufGraph <- load_graph("Buffalo_15380_invs.csv", NULL, NULL)

starts <- c(1975)#seq(1975, 2002)
window <- 6
Y <- list()
techHash <- new.env(hash=TRUE, parent=empty()) 
for(start_year in starts){
    g1 <- subset.on.edgeset(bufGraph,
                        which(E(bufGraph)$appyear >=start_year &
                            E(bufGraph)$appyear <=start_year+window),
                        TRUE)
    sapply((1:vcount(g1)) 
    adj <- get.adjacency(g1, binary=TRUE, type="lower")
    Y[[start_year]] <- adj[lower.tri(adj)]
    print(sum(Y[[start_year]]));
}

