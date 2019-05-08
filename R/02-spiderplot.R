#' Draw a haplotype bifurcation diagram
#'
#' @param hh an object of class `haplohh` (from rehh package)
#' @param focal genomic position of the focal site
#' @param left extend haplotypes this many (variant) sites to the left of the focal site
#' @param right extend haplotypes this many (variant) sites to the right of the focal site
#' @param palette name of an RColorBrewer palette used for colors
#' @param reverse logical; if `TRUE`, put derived allele in upper panel
#' @param nucleotides logical; if `TRUE`, put nucleotide label at each node to show how haplotype core changes
#' @param pmapobj logical; if \(nucleotides) is set to `TRUE`, a (pmap) object from the .inp file is required
#' @param nucleotidetrim_left numeric; positions of nucleotide annotations (left truncates)
#' @param nucleotidetrim_right numeric; positions of nucleotide annotations (right truncates)
#' @param scale_text numeric; a scaling factor for the size of the nucleotide text at the nodes. Must be between 0-1.
#' @param relabel a named vector to use for re-labelling panels; default labels are 'Ancestral' and 'Derived'
#' @param ... other arguments passed through to rehh::bifurcation.diagram()
#'
#' @value bifurcation diagram drawn with ggplot2
#'
#' @details needs ggplot2 and dplyr loaded in your environment; requires rehh and RColorBrewer, but only via a namespace

spiderplot <- function(hh, focal, nucleotides = T, pmapobj = NULL, nucleotidetrim_left = NA, nucleotidetrim_right = NA, scale_text = 0.5,
                       left = 10, right = 10, max.haps = 10, palette = "RdBu", reverse = FALSE, relabel = NULL, ...) {

  if(nucleotides == T & is.null(pmapobj)){
    stop("For the nucleotides options to work, you must specify a pmapobj")
    }

  ## find index of site nearest the desired focal position
  ifocal <- findInterval(focal, hh@position)
  pfocal <- hh@position[ifocal]

  ## get line segments for bifurcation plot
  ancestral = .my.bifurcation.diagram(hh, mrk_foc = ifocal, all_foc = 1, nmrk_l = left, nmrk_r = right, limhapcount = max.haps, ... )
  derived   = .my.bifurcation.diagram(hh, mrk_foc = ifocal, all_foc = 2, nmrk_l = left, nmrk_r = right, limhapcount = max.haps, ... )
  rez <- bind_rows(ancestral$rez, derived$rez,
                   .id = "allele" )

  rez$allele <- factor(rez$allele, levels = c(1,2), labels = c("ancestral", "derived"))

  ## find position of 'origin' -- place in middle of plot where all the bifurcations start
  find_origin <- function(df) {
    o <- dplyr::arrange(subset(df, x0 == pfocal), desc(n))$y0[1]
    tibble(origin = o, pos = pfocal)
  }
  origin <- group_by(rez, allele) %>% do(find_origin(.))

  ## set thickness of line segments proportional to (binned) haplotype frequency
  rez <- group_by(rez, allele) %>%
    mutate(bin = factor(as.numeric(cut(n, 5)) + 5*as.numeric(allele == "derived"), levels = 1:10))


  ## section for haplotype nodes
  hapnodes <- bind_rows(ancestral$hapnodes, derived$hapnodes,
                        .id = "allele" )

  hapnodes$allele <- factor(hapnodes$allele, levels = c(1,2),
                            labels = c("ancestral", "derived"))

  ## for some reason, these floats are different despite coming from same source, coerce with round
  pmapobj <- pmapobj %>%
    dplyr::rename(X = pos) %>%
    dplyr::mutate(X = round(X, 5))

  hapnodes <- hapnodes %>%
    dplyr::mutate(X = round(X, 5))

  hapnodes <- hapnodes %>%
    dplyr::left_join(x = ., y = pmapobj, by = c("X")) %>%
    dplyr::mutate(
      allelebase = ifelse(allelepol == 1, anc, ifelse(allelepol == 2, der, NA)),
      n = as.numeric(n),
      ntxt_scale = as.numeric(n)*scale_text)

  if(any(is.na(hapnodes$allelebase))){
    stop("There was an issue merging the .inp file with the EHH calculations for the \n
             ancestral or derived haplotypes")
  }

  if(!is.na(nucleotidetrim_right)){
    hapnodes <- hapnodes %>%
      dplyr::filter(X <= nucleotidetrim_right)
  }

  if(!is.na(nucleotidetrim_left)){
    hapnodes <- hapnodes %>%
      dplyr::filter(X >= nucleotidetrim_left)
  }


  ## set colors
  base.cols <- RColorBrewer::brewer.pal(11, palette)
  if (reverse)
    base.cols <- rev(base.cols)
  cols <- base.cols[ c(7:11,5:1) ]

  ## set default panel labels if none specified
  if (is.null(relabel))
    relabel <- c("ancestral"="Ancestral", "derived"="Derived")

  ## draw plot
  plotObj <- ggplot(rez) +
    geom_segment(aes(x = x0, xend = x1, y = y0, yend = y1, colour = bin, lwd = n), lineend = "round") +
    geom_point(data = origin, aes(x = pos, y = origin), size = 4, pch = 21, fill = "white") +
    scale_x_continuous("position (cM)") +
    scale_size_continuous(range = c(1,4)) +
    scale_colour_manual(values = cols) +
    guides(size = FALSE, colour = FALSE) +
    facet_grid(allele ~ ., scale = "free_y", space = "free_y",
               labeller = labeller(allele = relabel)) +
    theme(axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
          )

  if(nucleotides == T){
    plotObj <- plotObj +
      geom_point(data = hapnodes, aes(x = X, y = Y, size = n), pch = 21, fill = "#f0f0f0") +
      geom_text(data = hapnodes, aes(x = X, y = Y, label = allelebase, size = ntxt_scale), fontface = "bold")
  }

  return(plotObj)


  }

#' Get coordinates for lines in bifurcation diagram
#' @details Taken almost verbatim from rehh::bifurcation.diagram().
#' Not sure exactly how this does what it does, but seems to work.

.my.bifurcation.diagram <- function(haplohh,
                                    mrk_foc,
                                    all_foc,
                                    nmrk_l=2,
                                    nmrk_r=2,
                                    limhapcount=10,
                                    refsize=0.1,
                                    linecol="blue",
                                    main_leg=NA,
                                    xlab_leg="Position"){

  if(!(rehh:::is.haplohh(haplohh))){stop("Data oject is not of valid haplohh object... (see data2haplohh() function)")}
  if(nmrk_l<0 | nmrk_r<0){stop("nmrk_l and nmrk_r must be positive or null")}
  if(limhapcount<1){stop("limhapcount must be >1")}

  if(is.na(main_leg)){
    if(all_foc==1){main_leg=paste(haplohh@snp.name[mrk_foc], " (Ancestral Allele)",sep="")}
    if(all_foc==2){main_leg=paste(haplohh@snp.name[mrk_foc], " (Derived Allele)",sep="")}
  }
  #Checks
  if(mrk_foc-nmrk_l<0){
    stop("Too much markers on the left")
  }
  if(mrk_foc+nmrk_r>haplohh@nsnp){
    stop("Too much markers on the right")
  }

  #Recup haplos on the Right
  if(nmrk_r>0){
    haplo_r=haplohh@haplo[haplohh@haplo[,mrk_foc]==all_foc,][,(mrk_foc+1):(mrk_foc+nmrk_r)]
    haplo_r=haplo_r[rowSums(haplo_r==0)==0,] # drop any missing
    if(nrow(haplo_r)<limhapcount){stop("Number of available haplotypes on the right lower than limhapcount")}
    haplo_der_r=matrix(0,0,4) ; det_haplo_r=matrix(0,1,4) #; haplo_anc=matrix(0,0,4)
    colnames(haplo_der_r)=c("CODE_HAP","NhaploDer","HAP_DER1","HAP_DER2") #; colnames(haplo_anc)=c("CODE_HAP","HAP_ANC","n","mrk")
    colnames(det_haplo_r)=c("CODE_HAP","HAPLO","mrk","n")
    det_haplo_r[1,]=as.matrix(c(1,"",0,nrow(haplo_r)))

    for(mrk in 1:nmrk_r){
      nhaplo_tot=nrow(det_haplo_r)
      if(mrk==1){
        tmp=haplo_r[,1]
      }else{
        tmp=apply(haplo_r[,1:mrk],1,paste,collapse="")
      }
      list_hap_0=matrix(det_haplo_r[as.numeric(det_haplo_r[,3])==mrk-1,],ncol=4) #haplo au mrk n-1
      for(i in 1:nrow(list_hap_0)){
        hap_der=list()
        hap_der[[1]]=paste(list_hap_0[i,2],"1",sep="") ; hap_der[[2]]=paste(list_hap_0[i,2],"2",sep="")
        #print(hap_der)
        nhap_der=c(sum(tmp==hap_der[[1]]),sum(tmp==hap_der[[2]]))
        tmp_nhap_der=sum(nhap_der>0) ; tmp_hap_der=rep(0,2) ; tmp_cnt=0
        for(all in 1:2){
          if(nhap_der[all]>0){
            nhaplo_tot=nhaplo_tot+1
            det_haplo_r=rbind(det_haplo_r,c(nhaplo_tot,hap_der[[all]],mrk,nhap_der[all])) #; haplo_anc=rbind(haplo_anc,c(nhaplo_tot,list_hap_0[i,1],nhap_der[all],mrk))
            tmp_cnt=tmp_cnt+1 ; tmp_hap_der[tmp_cnt]=nhaplo_tot
          }
        }
        haplo_der_r=rbind(haplo_der_r,c(list_hap_0[i,1],tmp_nhap_der,tmp_hap_der))
      }
    }

    rownames(haplo_der_r)=haplo_der_r[,1] # ; rownames(haplo_anc)=haplo_anc[,1]
    haplo_der_r=haplo_der_r[,-1] #; haplo_anc=haplo_anc[,-1]
    #calcul coordonnees des haplo
    coord_r=matrix(0,nrow(det_haplo_r),2) ; rownames(coord_r)=det_haplo_r[,1] ; colnames(coord_r)=c("X","Y")
    for(i in nmrk_r:0){
      tmp_haplo=det_haplo_r[as.numeric(det_haplo_r[,3])==i,1]
      coord_r[tmp_haplo,1]=haplohh@position[i+mrk_foc]
      if(i==nmrk_r){
        coord_r[tmp_haplo,2]=1:length(tmp_haplo)/2
      }else{
        for(hap in tmp_haplo){
          if(as.numeric(haplo_der_r[hap,1])==1){coord_r[hap,2]=coord_r[haplo_der_r[hap,2],2]}
          if(as.numeric(haplo_der_r[hap,1])==2){coord_r[hap,2]=mean(coord_r[haplo_der_r[hap,2:3],2])}
          coord_r[hap,1]=haplohh@position[i+mrk_foc]
        }
      }
    }
  }

  # nfb addition
  # find inflection points, where are SNPs that cause the haplotypes to break
  inflxnpnt_right = as.data.frame(coord_r) %>%
    dplyr::group_by(X) %>% # X here is the original position where the SNP is location
    dplyr::summarise(poscounts = length(X)) %>% # count number of times X is repeated
    dplyr::filter(poscounts > 1) %>% # one is start, drop
    dplyr::group_by(poscounts) %>% # now find the min for each of these counts. The first position is the place where the "SNP" inflection points happens (e.g. we have to split from one line to two line and then from two to three, etc)
    dplyr::summarise(inflxnpnt = min(X)) %>% # min since we are moving upstream
    dplyr::select(c("inflxnpnt")) %>%
    unlist(.)

  # now find the coordinates of those inflection points and extract the alleles
  node_coords_haps_right <- cbind.data.frame(coord_r, det_haplo_r, side = "right", stringsAsFactors = F) %>%
    dplyr::filter(X %in% inflxnpnt_right) %>%
    dplyr::mutate(allelepol = base::substring(HAPLO, mrk))


  #Recup haplos on the Left
  if(nmrk_l>0){
    haplo_l=haplohh@haplo[haplohh@haplo[,mrk_foc]==all_foc,][,(mrk_foc-1):(mrk_foc-nmrk_l)]
    haplo_l=haplo_l[rowSums(haplo_l==0)==0,]
    if(nrow(haplo_l)<limhapcount){stop("Number of available haplotypes on the left lower than limhapcount")}
    haplo_der_l=matrix(0,0,4) ; det_haplo_l=matrix(0,1,4) #; haplo_anc=matrix(0,0,4)
    colnames(haplo_der_l)=c("CODE_HAP","NhaploDer","HAP_DER1","HAP_DER2") #; colnames(haplo_anc)=c("CODE_HAP","HAP_ANC","n","mrk")
    colnames(det_haplo_l)=c("CODE_HAP","HAPLO","mrk","n")
    det_haplo_l[1,]=as.matrix(c(1,"",0,nrow(haplo_l)))

    for(mrk in 1:nmrk_l){
      nhaplo_tot=nrow(det_haplo_l)
      if(mrk==1){
        tmp=haplo_l[,1]
      }else{
        tmp=apply(haplo_l[,1:mrk],1,paste,collapse="")
      }
      list_hap_0=matrix(det_haplo_l[as.numeric(det_haplo_l[,3])==mrk-1,],ncol=4) #haplo au mrk n-1
      for(i in 1:nrow(list_hap_0)){
        hap_der=list()
        hap_der[[1]]=paste(list_hap_0[i,2],"1",sep="") ; hap_der[[2]]=paste(list_hap_0[i,2],"2",sep="")
        nhap_der=c(sum(tmp==hap_der[[1]]),sum(tmp==hap_der[[2]]))
        tmp_nhap_der=sum(nhap_der>0) ; tmp_hap_der=rep(0,2) ; tmp_cnt=0
        for(all in 1:2){
          if(nhap_der[all]>0){
            nhaplo_tot=nhaplo_tot+1
            det_haplo_l=rbind(det_haplo_l,c(nhaplo_tot,hap_der[[all]],mrk,nhap_der[all])) #; haplo_anc=rbind(haplo_anc,c(nhaplo_tot,list_hap_0[i,1],nhap_der[all],mrk))
            tmp_cnt=tmp_cnt+1 ; tmp_hap_der[tmp_cnt]=nhaplo_tot
          }
        }
        haplo_der_l=rbind(haplo_der_l,c(list_hap_0[i,1],tmp_nhap_der,tmp_hap_der))
      }
    }

    rownames(haplo_der_l)=haplo_der_l[,1] # ; rownames(haplo_anc)=haplo_anc[,1]
    haplo_der_l=haplo_der_l[,-1] #; haplo_anc=haplo_anc[,-1]
    #calcul coordonnees des haplo
    coord_l=matrix(0,nrow(det_haplo_l),2) ; rownames(coord_l)=det_haplo_l[,1] ; colnames(coord_l)=c("X","Y")
    for(i in nmrk_l:0){
      tmp_haplo=det_haplo_l[as.numeric(det_haplo_l[,3])==i,1]
      coord_l[tmp_haplo,1]=haplohh@position[mrk_foc-i]
      if(i==nmrk_l){
        coord_l[tmp_haplo,2]=1:length(tmp_haplo)/2
      }else{
        for(hap in tmp_haplo){
          if(as.numeric(haplo_der_l[hap,1])==1){coord_l[hap,2]=coord_l[haplo_der_l[hap,2],2]}
          if(as.numeric(haplo_der_l[hap,1])==2){coord_l[hap,2]=mean(coord_l[haplo_der_l[hap,2:3],2])}
          coord_l[hap,1]=haplohh@position[mrk_foc-i]
        }
      }
    }
  }

  # nfb addition
  coord_l[,2]=coord_l[,2] + coord_r[1,2] - coord_l[1,2] # from below to get left side to slope down always

  # find inflection points, where are SNPs that cause the haplotypes to break for LEFT
  inflxnpnt_left = as.data.frame(coord_l) %>%
    dplyr::group_by(X) %>% # X here is the original position where the SNP is location
    dplyr::summarise(poscounts = length(X)) %>% # count number of times X is repeated
    dplyr::filter(poscounts > 1) %>% # one is start, drop
    dplyr::group_by(poscounts) %>% # now find the min for each of these counts. The first position is the place where the "SNP" inflection points happens (e.g. we have to split from one line to two line and then from two to three, etc)
    dplyr::summarise(inflxnpnt = max(X)) %>% # max now since we are moving downstream
    dplyr::select(c("inflxnpnt")) %>%
    unlist(.)

  # now find the coordinates of those inflection points and extract the alleles
  node_coords_haps_left <- cbind.data.frame(coord_l, det_haplo_l, side = "left", stringsAsFactors = F) %>%
    dplyr::filter(X %in% inflxnpnt_left) %>%
    dplyr::mutate(allelepol = base::substring(HAPLO, mrk))

  node_coords_haps <- rbind.data.frame(node_coords_haps_left, node_coords_haps_right, stringsAsFactors = F)



  #PLOT
  if(nmrk_l>0 & nmrk_r>0){
    coord_l[,2]=coord_l[,2] + coord_r[1,2] - coord_l[1,2]
    dum_coord=rbind(coord_l,coord_r)
  }
  if(nmrk_l==0){dum_coord=coord_r}
  if(nmrk_r==0){dum_coord=coord_l}

  #plot(dum_coord,pch="",yaxt="n",bty="n",xlab=xlab_leg, main=main_leg,ylab="")
  #abline(v=haplohh@position[mrk_foc],lty=2,col=linecol)

  rez <- NULL
  if(nmrk_r>0){
    #lwd_adjust=as.numeric(det_haplo_r[,4])*refsize/(as.numeric(det_haplo_r[1,1]))
    lwd_adjust <- as.numeric(det_haplo_r[,4])
    names(lwd_adjust)=rownames(coord_r)
    for(i in (nmrk_r-1):0){
      tmp_haplo=det_haplo_r[as.numeric(det_haplo_r[,3])==i,1]
      for(hap in tmp_haplo){
        x0=coord_r[hap,1] ; y0=coord_r[hap,2]
        for(j in 1:as.numeric(haplo_der_r[hap,1])){
          tmp_lwd=lwd_adjust[haplo_der_r[hap,1+j]]
          x1=coord_r[haplo_der_r[hap,1+j],1] ; y1=coord_r[haplo_der_r[hap,1+j],2]
          #lines(c(x0,x1),c(y0,y1),lwd=tmp_lwd,lty=1,col=linecol)
          rez <- rbind(rez, data.frame(x0 = x0, x1 = x1, y0 = y0, y1 = y1, side = "right",
                                       n = tmp_lwd))
        }
      }
    }
  }

  if(nmrk_l>0){
    #lwd_adjust=as.numeric(det_haplo_l[,4])*refsize/(as.numeric(det_haplo_l[1,1]))
    lwd_adjust <- as.numeric(det_haplo_l[,4])
    names(lwd_adjust)=rownames(coord_l)
    for(i in (nmrk_l-1):0){
      tmp_haplo=det_haplo_l[as.numeric(det_haplo_l[,3])==i,1]
      for(hap in tmp_haplo){
        x0=coord_l[hap,1] ; y0=coord_l[hap,2]
        for(j in 1:as.numeric(haplo_der_l[hap,1])){
          tmp_lwd=lwd_adjust[haplo_der_l[hap,1+j]]
          x1=coord_l[haplo_der_l[hap,1+j],1] ; y1=coord_l[haplo_der_l[hap,1+j],2]
          rez <- rbind(rez, data.frame(x0 = x0, x1 = x1, y0 = y0, y1 = y1, side = "left",
                                       n = tmp_lwd))
        }
      }
    }
  }

  return(list(rez = rez,
              hapnodes = node_coords_haps)
  )

}

