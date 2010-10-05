#This file includes some functions which are useful when building panels.
#There are also commented lines of code which you can run in the intepreter
#to see how the functions work. They're commented out so you can load the functions
#quickly.

#You can load this file into an R interpreter workspace by dropping this file into
#R's working directory (or switching R's working directory to this file's directory)
#and calling "source("panel_funcs.R")".

library(igraph)

#load("final_graph_full_msa.RData")
#Loaded graph will be called "g".

#Use this function to slice a graph by edge indices.
subset.on.edgeset <- function(g, edge_indxs, keep_singletons=F){
    drop_idx <- which(is.na(match(seq(ecount(g)), edge_indxs)))
    #igraph functions that take edge indices indexed from 0, but R indexes from 1,
    #so subtract 1 from indices returned from calls to R functions like "which" and "match".
    g <- delete.edges(g, drop_idx-1)
    if(!keep_singletons){
        vert_idx <- unique(as.vector(get.edges(g, E(g))))
        g <- subgraph(g, vert_idx)
    }
    return(g)
}

######################################
#Example year-window panel analysis
######################################

#Create 3-year windows from 1975 to 2008
#window_starts <- seq(1975, 2006, by=3)

#Example function that grabs edge indices by their appyear property, then slices the graph
#to only include those edges and delete all resulting singletons.
#Write functions like these that can build a query from an element of a list, then feed
#these functions to "lapply" to create a panel of graphs in one line.
choose_window <- function(g, window_start){
    #Access edge attributes using E(g)$attribute syntax.

    #View all edge attributes using list.edge.attributes(g).

    #You could also do this using E(g)[appyear >= window_start & appyear < window_start+3]
    #but this is _much_ slower.
    window_edges <- which(E(g)$appyear >= window_start & E(g)$appyear < window_start+3)
    return(subset.on.edgeset(g, window_edges, F))
}

#Create year-window panel data
#g_list <- lapply(window_starts, function(x){ choose_window(g, x) })

#Tack on names so you can access these lists easily
#names(g_list) <- window_starts

#Take a look at your network
#E(g_list[["1975"]])$appyear[1:10]

#You can also lapply analyses across g_list. For example, add eigenvector centrality
#as a vertex attribute to each graph.
#lapply(g_list, function(g){ V(g)$evcent <- evcent(g)$vector; return(g) })

#####################################
#Example MSA panel analysis
#####################################

#Read in the ZIP -> MSA mapping file
t <- read.table("zip_cbsa_200403.txt", fill=T,
        colClasses=c("character", "character"),
        col.names=c("zip", "msa"))

zip_to_msa <- function(t, zip){
    return(t$msa[match(zip, t$zip)])
}

add_msa <- function(g, t){
    g<-set.edge.attribute(g, "h_msa", E(g), zip_to_msa(t,E(g)$h_zipcode))
    g<-set.edge.attribute(g, "t_msa", E(g), zip_to_msa(t,E(g)$t_zipcode))
    return(g)
}

#g <- add_msa(g, t)
#save(g, file="final_graph_full_msa.RData")

#Example function that grabs edges where both inventors are in the MSA of interest and
#slices the graph to include only inventors incident to those edges.
choose_msa <- function(g, msa){
    msa_edges <- which(E(g)$h_msa == msa & E(g)$t_msa == msa)
    return(subset.on.edgeset(g, msa_edges, F))
}

#Create an MSA panel
#Warning: this line as written takes a while to run. You might want to run this on 
#subsets of MSAs, then merge them together into the full list. To to this, create
#a list of MSAs, and substitute that list in for the 'unique(t$msa)' here.

#g_list <- lapply(unique(t$msa), function(x){ choose_msa(g, x) })
#names(g_list) <- unique(t$msa)
