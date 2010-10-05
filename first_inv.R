library(RMySQL)
library(igraph)

##These top functions are used to create a network from my particular
##MySQL database. Currently configured for the MARA study dataset.

close_out_MySQL <- function(drv){
	sapply(dbListConnections(drv), function(x){ sapply(dbListResults(x), dbClearResult); dbDisconnect(x);});
}

get_queries <- function(query_set, window=NULL){
	return(switch(query_set,
		list(vquery="select invnum_N as name, firstname, lastname, loc from inventors5c group by invnum_N",
		     equery="select * from inv_edges"),
        new_edges=list(vquery="select invnum_N as name, name as textname from inv_verts",
                 equery="select * from inv_edges_new"),
        edges_0002=list(vquery="select invnum_N as name, name as textname, loc, city, zipcode, assignee from verts_0002",
                 equery="select * from edges_0002"),
        edges_0305=list(vquery="select invnum_N as name, name as textname, loc, city, zipcode, assignee from verts_0305",
                 equery="select * from edges_0305"),
        edges_0608=list(vquery="select invnum_N as name, name as textname, loc, city, zipcode, assignee from verts_0608",
                 equery="select * from edges_0608"),
		preMara=list(vquery="select invnum_N as name, name as textname, loc, zipcode, assignee from premara_verts",
			     equery="select * from premara_edges;"),
        preMaraNoAdj=list(vquery="select invnum_N as name, name as textname, loc, zipcode, assignee from premara_verts_noadj",
                 equery="select * from premara_edges_noadj;"),
        preMaraUsOnly=list(vquery="select invnum_N as name, name as textname, loc, zipcode, assignee from premara_verts_usonly",
                 equery="select * from premara_edges_usonly;"),
        preMaraLoc = list(vquery="select invnum_N as name, name as textname, u_locs from premara_verts_loc",
                 equery="select * from premara_edges_noadj;"),
        postMaraNoAdj=list(vquery="select invnum_N as name, name as textname, loc, zipcode, assignee from postmara_verts_noadj",
                 equery="select * from postmara_edges_noadj;"),
        postpostMaraNoAdj=list(vquery="select invnum_N as name, name as textname, loc, zipcode, assignee from postpostmara_verts_noadj",
                 equery="select * from postpostmara_edges_noadj;"),
		MI=list(vquery="select invnum_N as name, firstname, lastname, loc from MI_dataset group by invnum_N",
			equery="select * from edges_MI_real"),
        window=list(vquery=sprintf("select * from window_verts where window=%d", window),
                    equery=sprintf("select * from window_edges where window=%d", window))));
}

calc_premara_n <- function(g, n_nums, listify=T){
    #function designed for g defined by create_network("preMaraLoc")
    if(listify)
        V(g)$u_locs_vec <- sapply(V(g)$u_locs, function(x){strsplit(x, '|', fixed=T)})

    con <- dbConnect(MySQL(), user="alex", db="inventors")
    on.exit(close_out_MySQL(MySQL()))
    invs <- dbGetQuery(con, "select distinct(invnum_N) as invnum_N from premara5")
    inv_idx <- match(invs$invnum_N, V(g)$name)-1
    print("Have inventors!")

    res <- invs 
    
    for(n in n_nums){
        hoods <- neighborhood(g, order=n, nodes=inv_idx)
        print("Have hoods!")
        loc_list <- V(g)$u_locs_vec
        i <- 0
        counts <- lapply(hoods, function(h) {i <<- i+1; if(!i%%500) cat('.'); get_counts(loc_list, h)}) 
        counts <- matrix(unlist(counts),nc=3,byrow=T)
        print("Have counts!")
        tmp <- invs
        tmp$coinvs <- counts[,1] 
        tmp$num_nenf <- counts[,2] 
        tmp$num_nc <- counts[,3] 

        res <- merge(res, tmp, by="invnum_N", suffixes=c("",sprintf("%d",n)))
    }

    return(res)
}

fix_graphml <- function(filename){
    g <- read.graph(filename, format="graphml")
    graphml <- xmlTreeParse(filename)
    graph_props <- graphml$doc$children$graphml$children$graph
    graph_props_new <- addAttributes(graph_props, parse.nodes=vcount(g), parse.edge=vcount(e))
    graphml$doc$children$graphml$children$graph <- graph_props_new
    root <- graphml
    saveXML(root, graphml)
}

    
get_counts <- function(loc_list, hood){
    nenf <- c('AK', 'CA', 'CT', 'MN', 'MT', 'ND', 'NV', 'OK', 'WA', 'WV', 'MI')
    coinvs <- length(hood)-1

    u_loc <- unique(unlist(loc_list[hood+1]))
    nenf <- length(intersect(nenf, u_loc))
    enf <- length(u_loc)-nenf

    return(c(coinvs, nenf, enf))
}
    
            
        
    

create_network <- function(query_set = "default", window=NULL){
	#e <- loadEdges()
	#print("loaded edges!")
	#v <- loadVertices()
	#v <- v[order(v$invnum_N),]
	#print("loaded vertices!")
#	e$t <- match(e$t, v)
#	print("matched tails!")
#	e$h <- match(e$h, v)
#	print("matched heads!")
	g <- graph.empty(directed=F)
	queries <- get_queries(query_set, window)

	con <- dbConnect(MySQL(), user="alex", db="inventors")
	on.exit(close_out_MySQL(MySQL()))

	print("beginning vertex query")
	v_time <- system.time(v <- dbSendQuery(con, queries$vquery))
	print(v_time)
	while(!dbHasCompleted(v)){
		tmp_v <- fetch(v, 1000000)
		gc()
		g <- add.vertices(g, dim(tmp_v)[1], attr=tmp_v)
	}
	rm(tmp_v)
	dbClearResult(v)
	gc()
	print("done with vertices!")
	print(vcount(g))
	
	print("beginning edge query")
	e_time <- system.time(e <- dbSendQuery(con, queries$equery))
	print(e_time)
	i <- 0
	while(!dbHasCompleted(e)){
		tmp_e <- fetch(e, 500000)
		gc()
		tmp_e$h <- match(as.character(tmp_e$h), V(g)$name)-1
		tmp_e$t <- match(as.character(tmp_e$t), V(g)$name)-1
		g <- add.edges(g, matrix(c(tmp_e$h, tmp_e$t), nr=2, byrow=T), attr = tmp_e[,c(-1,-2)])
		i <- i+1
		if(!(i %% 3)) print(ecount(g))
	}
	rm(tmp_e)
	dbClearResult(e)
	print("done with edges!")
	print(ecount(g))

	dbDisconnect(con)

	gc()
	return(g)
}

loadEdges <- function(){
	con <- dbConnect(MySQL(), user="alex", password="", dbname="inventors")
	e_res <- dbSendQuery(con, "select * from inv_edges;")
	success <- dbApply
	dbDisconnect(con)
	return(e)
}

loadVertices <- function(){
	con <- dbConnect(MySQL(), user="alex", password="", dbname="inventors")
	v <- dbGetQuery(con, "select invnum_N,loc from inventors group by invnum_N")
	dbDisconnect(con)
	return(v)
}

#Like subgraph, but takes an edgeset and returns a graph that only contains those edges
#and vertices incident to those edges.
subset.on.edgeset <- function(g,es){
	e <- get.edges(g, es)
	gc()
	vset <- V(g)[unique(c(e[,1],e[,2]))]
	return(subgraph(g,vset))
}

graph2df <- function(g){
    vatts <- list.vertex.attributes(g)
    vdf <- data.frame(sapply(list.vertex.attributes(g),
                      function(x){get.vertex.attribute(g, x, V(g))},
                      simplify=F))
    everts <- matrix(V(g)$name[get.edges(g, E(g))+1], nc=2)
    if(length(list.edge.attributes(g)))
        edf <- data.frame(h=everts[,1], t=everts[,2], sapply(list.edge.attributes(g),
                      function(x){get.edge.attribute(g, x, E(g))},
                      simplify=F))
    else
        edf <- data.frame(h=everts[,1], t=everts[,2])

    return( list(vertices=vdf, edges=edf) )
}

#Exports a graph object to tulip format.
graph2tulip <- function(g,file="output.tlp"){
	cat(file=file, "(tlp \"2.0\"\n");

	vertices <- seq(1,(vcount(g)))
	cat(file=file, ",****vertices\n", append=T)
	cat(file=file, "\t( nodes", vertices, ")\n", append=T)

	edge_mat <- get.edges(g, E(g))+1
	cat(file=file, ";****edges\n", append=T)
	cat(file=file, paste("\t(edge", E(g)+1, edge_mat[,1], edge_mat[,2], ")\n"), append=T)

	for(prop in list.vertex.attributes(g)){
		thetype <- if(typeof(get.edge.attribute(g,prop)) == "double") "double" else "string";
		cat(file=file, sprintf(";****%s vertex property\n", prop), append=T)
		cat(file=file, sprintf("\t(property 0 %s \"%s\"\n", thetype, prop), append=T)
		cat(file=file, paste("\t\t(node", vertices, sprintf("\"%s\"",get.vertex.attribute(g, prop)), ")\n"), append=T)
		cat(file=file, sprintf(") ;%s vertex property\n", prop),append=T)
	}

	for(prop in list.edge.attributes(g)){
		thetype <- if(typeof(get.edge.attribute(g, prop)) == "double") "double" else "string";
		cat(file=file, sprintf(";****%s edge property\n", prop), append=T)
		cat(file=file, sprintf("\t(property 0 %s \"%s\"\n", thetype, prop), append=T)
		cat(file=file, paste("\t\t(edge ", E(g)+1, sprintf(" \"%s\"", get.edge.attribute(g, prop)), ")\n"), append=T)
		cat(file=file, sprintf(") ;%s edge property\n", prop),append=T)
	}

	cat(file=file, ")", append=T)
}

#Returns the largest component of a graph as a graph object. Faster than decompose.graph
component <- function(g,num=1){
	k <- clusters(g)
	comp_num <- as.numeric(names(sort(tapply(k$membership, k$membership, length), decreasing=T)[num]))
	g1 <- subgraph(g, V(g)[(which(k$membership==comp_num)-1)])
	return(g1)
}

#This function doesn't do exactly what it's supposed to yet. But it does add a
#w property, which is the "weight" of a particular connection.
add.cascading.weights <- function(g){
	for(i in sort(unique(E(g)$appyear),decreasing=T))
		E(g)[appyear <= i]$w <- count.multiple(g, E(g)[appyear <= i])
	return(g)
}
