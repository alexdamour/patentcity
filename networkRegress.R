library(igraph)
source("networkUtils.R")
source("panel_funcs.R")

pois_loglike <- function(par, X, y){
    alpha <- X %*% par
    lambda <- exp(alpha)

    ll <- -lambda + k*log(lambda)
    return(sum(ll))
}

solve_fit <- function(g, model){
    X <- switch(model, null=data.frame(rep(1, ecount(g))), stop("specify a model"))
    y <- as.vector(get.adjacency(simplify(g), attr="weight"))
    par <- dim(X)[2]
    res <- optim(par, gr=NULL, X, y)
}

buildY <- function(g, start_year, window){
    g1 <- subset.on.edgeset(g,
                        which(E(g)$appyear >=start_year &
                            E(g)$appyear <=start_year+window),
                        TRUE)
    adj <- get.adjacency(g1, binary=TRUE, type="lower")
    return(adj[lower.tri(adj)])
}

buildX <- function(g, orig_tab, start_year){
    cbref <- build_cumulative_tech(orig_tab, start_year);
    indices <- combn((1:vcount(g))-1, 2)
    id_index <- rbind(V(g)$invnum_N[indices[1,]], V(g)$invnum_N[indices[2,]])
    X <- apply(id_index, 2, function(x){
                                length(intersect(unlist(cbref[[as.character(x[1])]])),
                                                 unlist(cbref[[as.character(x[2])]])) >0
                                })
    print(sum(X))
    return(X)
}
    

build_cumulative_tech <- function(orig_tab, start){
    o <- order(orig_tab$appyear, orig_tab$invnum_N, decreasing=T)
    orig_tab <- orig_tab[o,]
    orig_tab <- orig_tab[orig_tab$appyear < start,]
    class_bags <- tapply(orig_tab$classes, orig_tab$invnum_N,
                                function(x){ unique(unlist(
                                    sapply(strsplit(
                                               unlist(strsplit(as.character(x),'|',fixed=T)),
                                               '-',fixed=T
                                           ), function(x){x[1]}
                                    )
                                  ))
                                })
    return(class_bags)
}

sample_nonedges <- function(g, mult=1){
    edges <- get.edges(g, E(g))
    num_edges <- ecount(g)*mult
    print(paste("Sampling", num_edges, "non-edges."))
    num_verts <- vcount(g)
    num_to_go <- num_edges 
    new_edges <- data.frame(head=c(), tail=c())
    while(num_to_go){
        print(paste("Have", num_edges - num_to_go, "of", num_edges, "edges..."))
        heads <- sample((1:num_verts)-1, num_to_go, replace=T)
        tails <- sample((1:num_verts)-1, num_to_go, replace=T)
        new_edges_app <- data.frame(head=heads, tail=tails)
        new_edges_app <- new_edges_app[!are.connected(g, new_edges_app$heads, new_edges_app$tails),]
        new_edges <- rbind(new_edges, new_edges_app)
        num_to_go <- num_edges - dim(new_edges)[1]
    }
    print(paste("Sampled", dim(new_edges)[1], "of", num_edges, "edges."))
    return(new_edges)
}

component <- function(g,num=1){
	k <- clusters(g)
	comp_num <- as.numeric(names(sort(tapply(k$membership, k$membership, length), decreasing=T)[num]))
	g1 <- subgraph(g, V(g)[(which(k$membership==comp_num)-1)])
	return(g1)
}

inspect_components <- function(listname){
    load(paste(listname, "RData", sep='.'))
    graph_list <- get(listname)
    lapply(graph_list, function(g){ E(g)$weight <- count.multiple(g) })
    print(lapply(graph_list, function(g){ vcount(component(g, 1))/vcount(g) }))
}

fit_model <- function(mylist, filename, year_start){
    buf_list <- mylist
    df <- read.csv(filename, stringsAsFactors=F)
    g <- buf_list[[as.character(year_start)]]
    cbref <- build_cumulative(g, df, year_start)

    edges <- get.edges(g, E(g))
    nonedges <- sample_nonedges(g, 1)
    keymat <- cbind(V(g)$invnum_N[edges[,1]+1], V(g)$invnum_N[edges[,2]+1])
    keymat <- rbind(keymat, cbind(V(g)$invnum_N[nonedges[,1]+1], V(g)$invnum_N[nonedges[,2]+1]))
    X <- apply(keymat, 1, function(x){length(intersect(unlist(cbref[[as.character(x[1])]]),
                                                       unlist(cbref[[as.character(x[2])]])))})
    #X <- cbind(1, X)
    Y <- c(rep(1, dim(edges)[1]), rep(0, dim(nonedges)[1]))
    return(glm(Y~X,family=binomial(link="logit")))
}

fit_giant_model <- function(filename, start_range, window){
    df <- read.csv(filename, stringsAsFactors=F)
    g <- load_graph(filename, NULL, NULL)
    active_dyads <- get.edges(g, E(g))
    nonedges <- as.matrix(sample_nonedges(g, 1))
    all_dyads <- rbind(active_dyads, nonedges)
    keymat <- cbind(V(g)$invnum_N[all_dyads[,1]+1], V(g)$invnum_N[all_dyads[,2]+1])
    g1 <- subset.on.edgeset(g, which(E(g)$appyear >= 1975 & E(g)$appyear <= 1981), TRUE)
    print(dim(all_dyads))
    Y_last <- data.frame(Y=(as.numeric(apply(all_dyads, 1, function(x){are.connected(g1, x[1], x[2])}))));
    bigY <- data.frame(Y=NULL)
    bigX <- data.frame(ol=NULL, new=NULL, old=NULL)

    for(s in start_range){
        g1 <- subset.on.edgeset(g,
                        which(E(g)$appyear >=s &
                            E(g)$appyear <=s+window),
                        TRUE)

        cbref <- build_cumulative_tech(df, s)
        print(head(cbref))
        X_overlap <- apply(keymat, 1, function(x){length(intersect(unlist(cbref[[as.character(x[1])]]),
                                                           unlist(cbref[[as.character(x[2])]])))>0})
        Y <- data.frame(Y=as.numeric(apply(all_dyads, 1, function(x){are.connected(g1, x[1], x[2])})));
        print(length(Y))
        print("Y")
        X_new <- Y-Y_last > 0
        X_old <- Y-Y_last < 0
        X <- data.frame(ol=X_overlap, fresh=X_new, stale=X_old)
        print("X")
        print(dim(X))
        Year <- rep(as.character(s), dim(X)[1])

        bigY <- rbind(bigY, Y)
        bigX <- rbind(bigX, X)
        Y_last <- Y
    }
    print(dim(bigY))
    print(dim(bigX))
    dimnames(bigX)[[2]] <- c("ol", "fresh", "stale")
    dattab <- cbind(bigY, bigX)
    dattab$year <- as.factor(unlist(sapply(start_range, function(x){rep(x, dim(keymat)[1])})))
    print(head(dattab))
    return(glm(Y~ol+fresh+stale+year,data = dattab, family=binomial(link="logit")))
}

###########SCRIPT BEGIN##############
#window_starts <- seq(1975, 2002, by=1)
#create_save_panels(window_starts, 5, T)
#inspect_components("buffalo_list")
#inspect_components("boston_list")
#inspect_components("sival_list")

#remake_lists(seq(1975,2002), 6)
#load("buf_list.RData")
#df <- read.csv("Buffalo_15380_invs.csv", stringsAsFactors=F)
#cbs <- lapply(buf_list, function(g){build_cumulative(g, df, 1976)})
#
#g <- buf_list[[2]]
#cbref <- cbs[[2]]
#edges <- get.edges(g, E(g))
#nonedges <- sample_nonedges(g, 1)
#keymat <- cbind(V(g)$invnum_N[edges[,1]+1], V(g)$invnum_N[edges[,2]+1])
#keymat <- rbind(keymat, cbind(V(g)$invnum_N[nonedges[,1]+1], V(g)$invnum_N[nonedges[,2]+1]))
#X <- apply(keymat, 1, function(x){length(intersect(unlist(cbref[[as.character(x[1])]]),
#                                                   unlist(cbref[[as.character(x[2])]])))>0})
#X <- cbind(1, X)
#Y <- c(rep(1, dim(edges)[1]), rep(0, dim(nonedges)[1]))
#fit <- glm.fit(X,Y,family=binomial(link="logit"))
#load("buf_list.RData")
#fit_list <- lapply(seq(1976,2002), function(x){fit_model(buf_list, "Buffalo_15380_invs.csv", x)})

filename <- "Buffalo_15380_invs.csv"
ss <- calc_graph_stats(filename, seq(1975,2003), 7)

#fit <- fit_giant_model("Buffalo_15380_invs.csv", 1976:2002, 6)
