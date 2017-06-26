#' Print table of KEGG pathways
#' 
#' Links to corresponding entries on the KEGG website are added to the output.
#' 
#' @param x An object of class \code{topKEGG}
#' @param convert LOgical indicating whether gene IDs should be converted.
#' @param genes Vector of entrez IDs to highlight.
#' @param map An environment mapping entrez IDs to gene symbols.
#' @param A \code{data.frame} with kegg annotations to use. If this is \code{NULL}
#' annotations will be downloaded from KEGG.
#' @param tooltip Column to use for tooltips.
#' @param ... Additional arguments are passed to \code{pander}.
#' @export
printMD.topKEGG <- function(x, convert=FALSE, genes=NULL, map=NULL, annot=NULL, tooltip=c('up', 'down'), ...){
  pathways <- rownames(x)
  if(!is.null(genes)){
    if(is.null(annot)){
      species <- stringr::str_match(pathways[1], 'path:(\\S{3})\\d+')[2]
      annot <- getGeneKEGGLinks(species, convert = convert)
    }
    annot <- subset(annot, GeneID %in% genes)
    gene_annot <- lapply(pathways, function(path) subset(annot, PathwayID == path)$GeneID)
  }
  pathways <- stringr::str_replace(pathways, '^path:', 'http://www.kegg.jp/pathway/')
  if(!is.null(genes)){
    pathways <- paste(pathways, sapply(gene_annot, function(genes) paste(genes, collapse='+')), sep='+')
    if(!is.null(map)){
      gene_sym <- lapply(gene_annot, mget, envir=map, ifnotfound=NA)
    }
    for(tt in tooltip){
      if(tt %in% names(x)){
        x[[tt]] <- mapply(add_tooltip, text=x[[tt]], entrez=gene_annot, symbols=gene_sym, MoreArgs=list(type='gene-list'))
      }
    }
  }
  x$Pathway <- paste0('[', x$Pathway, '](', pathways, ')')
  x <- make_ref_links(x, 'Pathway')
  rownames(x) <- NULL
  pval <- which(stringr::str_detect(names(x), stringr::fixed('P.')))
  for(i in pval){
    x[[i]] <- printMD(x[[i]])
  }
  tab <- htmlTable::htmlTable(x, css.rgroup = "", css.rgroup.sep = "", css.tspanner = "",
                              css.tspanner.sep = "", css.total = "", css.table="",
                              css.cell = "", css.cgroup = "", 
                              css.class = "table table-striped table-bordered table-hover table-condensed")
  cat(tab)
}

#' Print table of enriched GO terms
#' 
#' Links to GO terms on AmiGO
#' @param x An object of class topGO.
#' @param de List of entrez IDs for differentially expressed genes.
#' @param go_map Environment mapping go terms to entrez IDS.
#' @param symbol_map Environment mapping entrez IDs to gene symbols.
#' @param ... Ignored
#' 
#' @export
printMD.topGO <- function(x, de, go_map, symbol_map, ...){
  terms <- paste0('http://amigo.geneontology.org/amigo/term/', rownames(x))
  x$Term <- paste0('[', x$Term, '](', terms, ')')
  x <- make_ref_links(x, 'Term')
  if(!missing(de)){
    for(grp in names(de)){
      if(grp %in% names(x)){
        ids <- mget(rownames(x), envir = go_map, ifnotfound = NA)
        ids <- mapply(intersect, ids, de[grp])
        sym <- lapply(ids, mget, envir=symbol_map, ifnotfound=NA)
        sym <- lapply(sym, unlist)
        x[[grp]] <- mapply(add_tooltip, text=x[[grp]], entrez=ids, symbols=sym, 
                           MoreArgs=list(type='gene-list'))
      }
    }
  }
  rownames(x) <- NULL
  pval <- which(stringr::str_detect(names(x), stringr::fixed('P.')))
  for(i in pval){
    x[[i]] <- printMD(x[[i]])
  }
  tab <- htmlTable::htmlTable(x, css.rgroup = "", css.rgroup.sep = "", css.tspanner = "",
                              css.tspanner.sep = "", css.total = "", css.table="",
                              css.cell = "", css.cgroup = "", 
                              css.class = "table table-striped table-bordered table-hover table-condensed")
  cat(tab)
}

#' @export
printMD.camera <- function(x, de_sets, map, de, expr, set_ids, groups, suffix, ...){
  x$pmid <- link_pubmed(x$pmid)
  if(!missing(de_sets)){
    de_sets <- de_sets[x$set_id]
    sym <- lapply(de_sets, mget, envir=map, ifnotfound=NA)
    sym <- lapply(sym, unlist)
    x[['DE genes']] <- mapply(add_tooltip, text=x[['DE genes']], entrez=de_sets,
                              symbols=sym, MoreArgs=list(type='gene-list'))
  }
  if(!missing(set_ids)){
    set_entrez <- lapply(set_ids[x$set_id], unlist)
    set_sym <- lapply(set_entrez, mget, envir=map, ifnotfound=NA)
    x[['NGenes']] <- mapply(add_tooltip, text=x[['NGenes']], entrez=set_entrez,
                            symbols=set_sym, MoreArgs=list(type='gene-list'))
  }
  if(!missing(de)){
    ## create and link enrichment plots
    if(missing(groups)){
      groups <- sapply(strsplit(colnames(expr), '_'), '[[', 1)
    }
    for(i in 1:nrow(x)){
      base_name <- x$set_id[i]
      if(!missing(suffix)){
        base_name <- paste(base_name, suffix, sep='_')
      }
      file_name <- file.path(knitr::opts_chunk$get('fig.path'), 
                             paste0(base_name, '.png'))
      if(!file.exists(file_name)){
        plot <- geneset_plot(direction=tolower(x[['Direction']][i]),
                             set=set_ids[[x$set_id[i]]],
                             expr=expr, de=de, groups=groups,
                             title=paste('Enrichment for gene set', x$set_id[i]),
                             ...)
        if(!dir.exists(knitr::opts_chunk$get('fig.path'))){
          dir.create(knitr::opts_chunk$get('fig.path'), recursive=TRUE)
        }
        png(file_name, width=knitr::opts_chunk$get('fig.width'), 
            height=knitr::opts_chunk$get('fig.height'), units='in', res=knitr::opts_chunk$get('dpi'))
        print(plot)
        dev.off()
      }
      x[['description']][i] <- add_tooltip(text=x[['description']][i], image=file_name,
                                           type='geneset-plot')
    }
  }
  
  x$FDR <- printMD(x$FDR)
  x$PValue <- printMD(x$PValue)
  x <- make_ref_links(x, which(names(x) == 'pmid'))
  tab <- htmlTable::htmlTable(x, css.rgroup = "", css.rgroup.sep = "", css.tspanner = "",
                              css.tspanner.sep = "", css.total = "", css.table="",
                              css.cell = "", css.cgroup = "", 
                              css.class = "table table-striped table-bordered table-hover table-condensed")
  cat(tab)
}

#' Print table of differential expression results
#' @param x An object of class \code{geneList}
#' @param ... Ignored
#' @export
printMD.geneList <- function(x, ...){
  if('AveExpr' %in% names(x)){
    last_col <- max(which(stringr::str_detect(names(x), '_confGroup$|P\\.Val|_logFC$')))
    if(last_col < ncol(x)){
      expr <- x[,(last_col + 1):ncol(x)]
      group <- stringr::str_match(names(expr), '(.+)[_.][^_.]+')[,2]
      x$AveExpr <- expr_tooltip(expr=expr, group=group, gene=x$symbol, text=printMD(x$AveExpr))
      x <- x[, -((last_col + 1):ncol(x))]
    }
  }
  x$symbol <- link_gene(x$entrez, x$symbol)
  x <- select(x, -entrez)
  for(i in which(sapply(x, is.numeric))){
    x[[i]] <- printMD(x[[i]])
  }
  names(x) <- stringr::str_replace(names(x), '_', ' ')
  tab <- htmlTable::htmlTable(x, css.rgroup = "", css.rgroup.sep = "", css.tspanner = "",
                              css.tspanner.sep = "", css.total = "", css.table="",
                              css.cell = "", css.cgroup = "", 
                              css.class = "table table-striped table-bordered table-hover table-condensed")
  cat(tab)
}

#' Create link to PubMed entry from PMID
#' @param pmid A vector of PMIDs
#' @param format Format of links to be created.
#' @export
link_pubmed <- function(pmid, format=c('markdown', 'html')){
  format <- match.arg(format)
  link <- paste0('https://www.ncbi.nlm.nih.gov/pubmed/', pmid)
  if(format == 'markdown'){
    link <- paste0('[', pmid, '](', link,')')
  } else{
    link <- paste0('<a href=', link,  '>', pmid, '</a>')
  }
  ifelse(pmid == '', '', link)
}

#'@export
link_gene <- function(entrez, text, format=c('markdown', 'html')){
  format <- match.arg(format)
  if(missing(text)){
    text <- entrez
  }
  link <- 'https://www.ncbi.nlm.nih.gov/gene?cmd=Retrieve&dopt=full_report&list_uids='
  if(format == 'markdown'){
    paste0('[', text, '](', link, entrez, ')')
  } else {
    paste0('<a href=', link, entrez, '>', text, '</a>')
  }
}

expr_tooltip <- function(expr, group, gene, text){
  if(missing(group)){
    group <- sapply(strsplit(names(expr), '_'), '[[', 1)
  }
  ans <- character(nrow(expr))
  for(i in 1:nrow(expr)){
    file_name <- file.path(knitr::opts_chunk$get('fig.path'), 
                           paste0(gene[i], '.png'))
    if(!file.exists(file_name)){
      plot <- expression_plot(unlist(expr[i, ]), group, main=gene[i])
      if(!dir.exists(knitr::opts_chunk$get('fig.path'))){
        dir.create(knitr::opts_chunk$get('fig.path'), recursive=TRUE)
      }
      png(file_name, width=knitr::opts_chunk$get('fig.width'), 
          height=knitr::opts_chunk$get('fig.height'), units='in', res=knitr::opts_chunk$get('dpi'))
      print(plot)
      dev.off()
    }
    ans[i] <- add_tooltip(text=text[i], image=file_name, type='geneset-plot')
  }
  ans
}

#' Augment a string with a tooltip
#' @export
add_tooltip <- function(text, tooltip, type='tooltip', ...){
  if(missing(tooltip)){
    content <- ''
  } else {
    content <- paste0('data-content=', paste(tooltip, collapse=', '))
  }
  opts <- list(...)
  if(length(opts)){
    opts <- lapply(opts, function(x) if(length(x) > 1) paste0('"', paste(x, collapse=','), '"') else x)
    names(opts) <- paste0('data-', names(opts))
    vals <- unlist(opts)
    vals <- ifelse(vals == '', names(opts), vals)
    opts <- paste0(names(opts), '=', vals)
    opts <- paste(opts, collapse=' ')
  } else {
    opts <- ''
  }
  paste0('<a data-toggle=\"', type, '\" href="#" title="" ', content, ' ', 
         opts, '>', text, '</a>') 
}