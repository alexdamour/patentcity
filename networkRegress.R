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

fit_logit_model <- function(df, graph_windows){
    bigY <- data.frame(Y=NULL)
    bigX <- data.frame(ol=NULL, new=NULL, old=NULL)

    offsets <- 0

    years <- names(graph_windows)
    g_last <- graph_windows[[years[1]]]
    edge_last_hash <- build_edge_env(g_last)
    e_last <- standard_edges(g_last)

    for(year in years[-1]){
        print(year)
        g <- graph_windows[[year]]
        e_cur <- standard_edges(g)
        edge_cur_hash <- build_edge_env(g)

        active_dyads <- get.edges(g, E(g))

        stale <- unlist(sapply(setdiff(e_last, e_cur), function(e) edge_last_hash[[e]]))
        just_expired <- get.edges(g_last, stale-1)
        print("built stale, expired")

        keymat_last <- cbind(V(g_last)$invnum_N[just_expired[,1]+1],
                     V(g_last)$invnum_N[just_expired[,2]+1])
        nonedges <- as.matrix(sample_nonedges(add_by_invnum(g, keymat_last), 4))
        cur_dyads <- rbind(active_dyads, nonedges)

        #Observations
        Y <- c(rep(1, dim(active_dyads)[1]), rep(0, dim(nonedges)[1]+dim(just_expired)[1]))
        print("Y")

        #Tech class overlap
        keymat <- cbind(V(g)$invnum_N[cur_dyads[,1]+1], V(g)$invnum_N[cur_dyads[,2]+1])

        keymat <- rbind(keymat, keymat_last)
        cbref <- build_cumulative_tech(df, as.numeric(year))
        #print(head(cbref))
        X_overlap <- apply(keymat, 1, function(x){
                length(intersect(unlist(cbref[[as.character(x[1])]]),
                                 unlist(cbref[[as.character(x[2])]])))>0
        })

        Y_last <- numeric(length(Y))
        print(head(setdiff(e_cur, e_last)))
        from_last <- sapply(intersect(e_cur, e_last), function(e) edge_cur_hash[[e]])
        print(head(from_last))
        Y_last[from_last] <- 1
        Y_last[-(1:dim(cur_dyads)[1])] <- 1

        Year <- rep(year, length(Y))

        X <- data.frame(ol=X_overlap, Y_last=Y_last, year=Year)
        print("X")
        print(dim(X))
        offsets <- c(offsets, length(Y))

        Y <- as.matrix(Y, nc=1)
        bigY <- rbind(bigY, Y)
        bigX <- rbind(bigX, X)

        g_last <- g
        e_last <- e_cur
        edge_last_hash <- edge_cur_hash
    }
    print(dim(bigY))
    print(dim(bigX))
    #dimnames(bigX)[[2]] <- c("ol", "fresh", "stale", "year")
    dattab <- cbind(bigY, bigX)
    return(list(dat=dattab, offsets=cumsum(offsets)))
    #print(head(dattab))
    #return(glm(Y~ol+fresh+stale+year,data = dattab, family=binomial(link="logit")))
}

rank_edges <- function(g_list, fit, offsets){
    y <- sort(names(g_list))[-1]
    for(i in 1:length(y)){
        g <- g_list[[y[i]]]
        E(g)$rank <- fit$fitted.values[offsets[i]+(1:ecount(g))]
        print(all(fit$y[offsets[i]+(1:ecount(g))]==1))
        g_list[[y[i]]] <- g
    }
    return(g_list)
}

component_percentages <- function(ranked_g_list, threshold){
    y <- sort(names(ranked_g_list))[-1]
    pct <- c()
    compsize <- c()
    real_compsize <- c()
    for(i in 1:length(y)){
        g0 <- ranked_g_list[[y[i]]]

        g_real <- component(g0, 1)
        real_compsize <- c(real_compsize, vcount(g_real))
        pct <- c(pct, length(E(g_real)$rank[E(g_real)$rank >= threshold])/ecount(g_real))

        g_fake <- delete.edges(g0, E(g0)[rank < threshold])
        compsize <- c(compsize, vcount(component(g_fake, 1)))
    }

    res <- data.frame(pct=pct, compsize=compsize, real_compsize=real_compsize)
    dimnames(res)[[1]] <- y
    return(res)
}

postprocess <- function(g_list, fit, offsets, filestub){
    gl <- rank_edges(g_list, fit, offsets)
    res50 <- component_percentages(gl, 0.5)
    res95 <- component_percentages(gl, 0.95)

    pdf(sprintf("%s_comp_vcount.pdf", filestub))
    plot(as.numeric(dimnames(res50)[[1]]), res50$compsize, 
             type='b', lwd=2, col="red",
             main = "Giant Component Size Test",
             ylab = "Giant Component Size in Vertices",
             xlab="Year", ylim=c(0, 300))
    lines(as.numeric(dimnames(res95)[[1]]), res95$compsize,
             type='b', lwd=2, col="blue")
    lines(as.numeric(dimnames(res2)[[1]]), res2$real_compsize,
             type='b', lwd=2)
    legend("topleft", legend = c("Real", "50%", "95%"), col=c("black", "red", "blue"),
            lwd=c(2, 2, 2))
    dev.off()

    pdf(sprintf("%s_comp_edge_pct.pdf", filestub))
    plot(as.numeric(dimnames(res50)[[1]]), res50$pct,
             type='b', lwd=2, col="red",
             main = "Giant Component Edge Acceptance Rate",
             ylab = "Percentage of Giant Component Edges Kept",
             xlab="Year", ylim=c(0, 1))
    lines(as.numeric(dimnames(res95)[[1]]), res95$pct,
             type='b', lwd=2, col="blue")
    legend("bottomleft", legend = c("50%", "95%"), col=c("red", "blue"),
            lwd=c(2, 2))
    dev.off()
}
        

###########SCRIPT BEGIN##############
filename <- "Buffalo_15380_invs.csv"
df <- read.csv(filename, stringsAsFactors=FALSE, header=TRUE)
print("Loaded df...")

ss <- calc_graph_stats(filename, seq(1989,2002), 7)
#print(ss$stat_df)
buffalo_list <- ss$graph_list
save(buffalo_list, file="buffalo_list.RData")
print("Created graph windows...")

#load("buffalo_list.RData")

datmat <- fit_logit_model(df, buffalo_list)
print("Created regression matrix...")
fit <- glm(V1 ~ ol + Y_last, data=datmat$dat, family=binomial(link="logit"))
print("Fitted model...")
postprocess(buffalo_list, fit, datmat$offsets, "gt1980")
print("Postprocessed.")
