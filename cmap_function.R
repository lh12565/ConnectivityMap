#function
a =function(V,n){
    V %<>% as.matrix
    t = nrow(V)
    if(nrow(V)==0){
        return(0)
    }
    ((1:t)/t - V/n) %>% apply(2,max)
}


b = function(V,n){
    V %<>% as.matrix
    t = nrow(V)
    if(nrow(V)==0){
        return(0)
    }
    (V/n - (1:t - 1)/t) %>% as.matrix%>% apply(2,max)
}

ks = function(a,b){
    (a>b)*a + (b>a)*(-b)
}

ksCalc = function(V,n){
    a = a(V,n)
    b = b(V,n)
    return(ks(a,b))
}

score = function(kUp,kDown){
    (sign(kUp)!=sign(kDown))*(kUp-kDown)
}

scoreCalc = function(Vup,Vdown,n){
    aUp = a(Vup,n)
    bUp = b(Vup,n)
    kUp = ks(aUp,bUp)
    
    aDown = a(Vdown,n)
    bDown = b(Vdown,n)
    kDown = ks(aDown,bDown)
    score = score(kUp,kDown)
    return(data.frame(kUp = kUp, kDown = kDown,score = score,instance= colnames(Vup)))
}


connectivityMapEnrichment = function(upTags,downTags,rankMatrix,instances,pAdjustMethod = 'fdr' ,d=100000){
    n = rankMatrix %>% nrow 
    Vup = rankMatrix[upTags,] %>% apply(2,sort)
    Vdown = rankMatrix[downTags,] %>% apply(2,sort) 
    scores = scoreCalc(Vup,Vdown,n)
    
    # "p to be max( si) and q to be min( si) across all instances in the collection c"
    p = max(scores$score)
    q = min(scores$score)
    
    scores %<>% dplyr::mutate(ConScore = (score>0)*(score/p) + (score<0)*(-score/q))
    
    # "The Kolmogorov-Smirnov statistic is computed for the set of t instances
    # in the list of all n instances in a result ordered in descending order of
    # connectivity score and up score (see how connectivity score is
    # calculated), giving an enrichment score ks0."
    scores = scores %>% dplyr::arrange(desc(ConScore),desc(kUp))
	scores$cmap_name<-instances[match(scores$instance,rownames(instances)),"cmap_name"]
    
    chems = instances$cmap_name %>% unique

    confidence = chems %>% sapply(function(chem){
        # browser()
        # print(chem)
        chemInstances = rownames(instances)[instances$cmap_name %in% chem]
        
        relevantInstances = scores %>% dplyr::filter(instance %in% chemInstances)
        
        relevantUpCount = (relevantInstances$score>0) %>% sum
        relevantDownCount = (relevantInstances$score<0) %>% sum
        relevantMehCount = (relevantInstances$score==0) %>% sum
        nonNull = (relevantUpCount>=relevantDownCount)*relevantUpCount + 
            (relevantUpCount<relevantDownCount)*(relevantDownCount)
        nonNull = nonNull/nrow(relevantInstances)  #Proportion of consistent directions within the same drug
        
        V = match(chemInstances,scores$instance) %>% sort
        ks0 = ksCalc(V,nrow(instances))
        
        Vrandoms = memoRandomV(length(chemInstances),d)
        ksPerm = memoKsCalc(Vrandoms,nrow(instances))
        
        q = sum(abs(ksPerm) >= abs(ks0))
        p = q/d
        
        return(c(enrichment = ks0,
                 p = p,
                 nonNull = nonNull,
                 instanceCount =  length(chemInstances)))
    }) %>% t
    
    confidence %<>% as.data.frame
    confidence$FDR =  p.adjust(confidence$p,method = pAdjustMethod)
    confidence = confidence[c('enrichment','p','FDR','nonNull','instanceCount')]
    
    return(list(instanceScores = as.data.frame(scores), chemScores = confidence))
}


score_function<-function(degs,log2FoldChange=0,deg_significance='pvalue',score_significance='p'){
	if (!deg_significance %in% c("pvalue", "padj")) {
		stop("deg_significance must be 'pvalue' or 'padj'")
	}
	if (!score_significance %in% c("p", "FDR")) {
		stop("score_significance must be 'p' or 'FDR'")
	}

	degs_f<-degs[ degs$log2FoldChange>log2FoldChange&degs[[deg_significance]]<0.05,]
	upTags<-degs_f %>%dplyr::arrange_(.dots = c(glue('desc(log2FoldChange)'))) %>%
		dplyr::select(gene_name) %>% unlist
		#%>%mouse2human %>% {.$humanGene} %>% unique
	upTags<-unique(gpl96[ gpl96$GeneSymbols%in%upTags,"Probe"])

	degs_f<-degs[ degs$log2FoldChange<log2FoldChange&degs[[deg_significance]]<0.05,]
	downTags<-degs_f %>%dplyr::arrange_(.dots = c(glue('log2FoldChange'))) %>%
		dplyr::select(gene_name) %>% unlist
		#%>%mouse2human %>% {.$humanGene} %>% unique
	downTags<-unique(gpl96[ gpl96$GeneSymbols%in%downTags,"Probe"])

	#connectivity score
	out = connectivityMapEnrichment(upTags,downTags, rankMatrix,instances,d = 100000)

	#specificity score
	out$chemScores$specificity = 1:nrow(out$chemScores) %>%  sapply(function(i){
		chem = rownames(out$chemScores)[i]
		enrichments = MSigDB_enrich %>% sapply(function(x){
			x[chem,'enrichment']
		})
		sign = sign(out$chemScores[chem,'enrichment'])
		if(sign==-1){
			enrichments = enrichments[enrichments<=0]
		} else if(sign ==1){
			enrichments = enrichments[enrichments>=0]
			
		}
		sum(abs(out$chemScores[chem,'enrichment']) < abs(enrichments))/length(enrichments)
	})
	
	#reliability score
	out$chemScores$reliable = out$chemScores$nonNull>0.5 & out$chemScores$instanceCount > 1 & out$chemScores[[score_significance]]<0.05
	out$hits =  out$chemScores[out$chemScores$reliable ==TRUE & out$chemScores[[score_significance]] <0.05,] %>% {.[order(.$p),]} %>% rownames
	return(out)
}


