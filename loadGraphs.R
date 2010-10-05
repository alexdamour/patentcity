library(igraph)
source("panel_funcs.R")

load_graph <- function(filename, start_year, window){
    g <- graph.empty(directed=F)
    df <- read.csv(filename, stringsAsFactors=F)
    drop_cols <- c("appdate", "gdate", "claims", "classes", "scls_cnt", "scls_1", "scls_pair_1", "back5", "back5_sc", "back_sc", "forw5", "forw5_sc", "forw", "forw_sc", "X")
    for(d in drop_cols)
        df[[d]] <- NULL
    if(!is.null(start_year) & !is.null(window))
        df <- df[df$appyear >=start_year & df$appyear <= (start_year + window),]
    verts <- aggregate(df, by=list(df$invnum_N), min)
    print(length(unique(df$invnum_N)))
    
    g <- add.vertices(g, dim(verts)[1], attr=data.frame(invnum_N=verts$invnum_N))
    print(list.vertex.attributes(g))
    edges <- merge(df, df, by=c("patent", "pat_type", "appyear", "assignee", "numasg"), suffixes=c("_h", "_t"))
    edges <- edges[edges$invnum_N_t != edges$invnum_N_h,]
    #return(list(edges, g))

    edges$h <- match(as.character(edges$invnum_N_h), V(g)$invnum_N)-1
    edges$t <- match(as.character(edges$invnum_N_t), V(g)$invnum_N)-1
    g <- add.edges(g, matrix(c(edges$h, edges$t), nr=2, byrow=T), attr=edges)
    g <- subset.on.edgeset(g, (E(g)+1), F)
 
    return(g)
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

