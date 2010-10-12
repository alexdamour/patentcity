library(igraph)
source("panel_funcs.R")

#Load a graph from a flat file, subsetting for a particular year window.
#Very important for picking up unconnected inventors.
load_graph <- function(df, start_year, window){
    g <- graph.empty(directed=F)
    #df <- read.csv(filename, stringsAsFactors=F)
    drop_cols <- c("appdate", "gdate", "claims", "classes", "scls_cnt", "scls_1", "scls_pair_1", "back5", "back5_sc", "back_sc", "forw5", "forw5_sc", "forw", "forw_sc", "X")
    for(d in drop_cols)
        df[[d]] <- NULL
    if(!is.null(start_year) & !is.null(window))
        df <- df[df$appyear >=start_year & df$appyear < (start_year + window),]
    verts <- aggregate(df, by=list(df$invnum_N), min)
    print(length(unique(df$invnum_N)))
    print(dim(verts))
    
    g <- add.vertices(g, dim(verts)[1], attr=data.frame(invnum_N=verts$invnum_N))
    print(list.vertex.attributes(g))
    edges <- merge(df, df, by=c("patent", "pat_type", "appyear", "assignee", "numasg"), suffixes=c("_h", "_t"))
    edges <- edges[edges$invnum_N_t != edges$invnum_N_h,]
    #return(list(edges, g))

    edges$h <- match(as.character(edges$invnum_N_h), V(g)$invnum_N)-1
    edges$t <- match(as.character(edges$invnum_N_t), V(g)$invnum_N)-1
    g <- add.edges(g, matrix(c(edges$h, edges$t), nr=2, byrow=T), attr=edges)
    #g <- subset.on.edgeset(g, (E(g)+1), F)
 
    return(g)
}

#Build an edge hashlist
build_edge_env <- function(g){
    e_env <- new.env(hash=T)
    keys <- standard_edges(g)
    sapply(1:ecount(g), function(i){ e_env[[keys[i]]] <<- i })
    e_env
}

#Standardize the edge representation to find it in the master hashlist
standard_edges <- function(g){
    g <- simplify(g)
    ge <- get.edges(g, E(g))
    invmap <- V(g)$invnum_N
    ge <- apply(ge, 1, function(x){ sprintf("%s-%s",
                                            invmap[min(x)+1],
                                            invmap[max(x)+1]) })
}

#Calculate "sufficient" statistcs for all of the window subgraphs in the
#date range
calc_graph_stats <- function(filename, start_range, window){
    df <- read.csv(filename, stringsAsFactors=F, header=T)
    g <- simplify(load_graph(df, NULL, NULL))
    active_dyads <- get.edges(g, E(g))
    rm(g)

    g0 <- load_graph(df, start_range[1], window)
    #g0.tmp <- subset.on.edgeset(g0, E(g0)+1, FALSE)
    oldEdges <- standard_edges(g0)
    giant <- vcount(component(g0, 1))/vcount(g0)
    stat_df <- data.frame(s=start_range[1], addE=NA, delE=NA, addN=NA, delN=NA,
                            edges=ecount(simplify(g0)), verts=vcount(g0), giant=giant)
    #g0 <- g0.tmp

    for(s in start_range[-1]){
        g1 <- load_graph(df, s, window) 

        #newEdges <- as.numeric(apply(active_dyads, 1, function(x){
        #            are.connected(g1, x[1], x[2])}))
        newEdges <- standard_edges(g1)
        addE <- length(setdiff(newEdges, oldEdges))
        delE <- length(setdiff(oldEdges, newEdges))
        oldEdges <- newEdges

#g0 <- subset.on.edgeset(g0, E(g0)+1, FALSE)
        #g1 <- subset.on.edgeset(g1, E(g1)+1, FALSE)
        addN <- length(setdiff(V(g1)$invnum_N, V(g0)$invnum_N))
        delN <- length(setdiff(V(g0)$invnum_N, V(g1)$invnum_N))

        giant <- vcount(component(g1,1))/vcount(g1)

        stat_df <- rbind(stat_df, data.frame(s=s, addE=addE, delE=delE, addN=addN, delN=delN,
                                        edges=ecount(simplify(g1)), verts=vcount(g1), giant=giant))
        g0 <- g1
    }
    stat_df
}

load_windows <- function(filename, starts, window){
    return(lapply(starts, function(x){ load_graph(filename, x, window) }))
}

#Returns the largest component of a graph as a graph object. Faster than decompose.graph
component <- function(g,num=1){
	k <- clusters(g)
	comp_num <- as.numeric(names(sort(tapply(k$membership, k$membership, length), decreasing=T)[num]))
	g1 <- subgraph(g, V(g)[(which(k$membership==comp_num)-1)])
	return(g1)
}

remake_lists <- function(starts, window){
    buf_list <<- load_windows("Buffalo_15380_invs.csv", starts, window) 
    names(buf_list) <- starts
    save(buf_list, file="buf_list.RData")

#    sival_list <<- load_windows("SiVal_41940_invs.csv", starts, window)
#    names(sival_list) <- starts
#    save(sival_list, file="sival_list.RData")

#    boston_list <<- load_windows("Boston_14460_invs.csv", starts, window)
#    names(boston_list) <- starts
#    save(boston_list, file="boston_list.RData")
}

