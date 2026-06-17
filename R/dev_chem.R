#' @description PAVE from MSdev.
#' @describeIn PVAE initialize from MSdev
#' @title PAVE analysis
#' @param object MSdev object
#'
#' @returns PAVE
#' @export
#'
#' @examples
#' get_PAVE_from_MSdev(msdev)
get_PAVE_from_MSdev <- function(object){

  object@sampleInfo <-
    object@sampleInfo%>%
    dplyr::mutate(
      isotope_tracer = paste0(ifelse(grepl("13C",sample.name),"[13]C","")),
      isotope_tracer = paste0(isotope_tracer,ifelse(grepl("15N",sample.name),"[15]N","")),
      .after = sample.name
    )
  message_with_time("Please check the column ",
                    crayon::red("isotope_tracer"),
                    " to confirm tracer")

  object <- MSdev_checkSampleInfo(object)
  object <- MSdev_update_xcms_pdata(object)
  return(object)


}


PAVE_get_atom_count <- function(object, BPPARAM = SnowParam(workers = 6,progressbar = T)){


  polarity.index <- c("0" = "Negative",
                      "1"="Positive")
  for (i.pol in 0:1) {


    time.start <- Sys.time()

    ### get xcms
    {

      polarity.tag <- paste0(polarity.index[as.character(i.pol)],"MS1")
      xcms.xcms <- object@xcmsData[[polarity.tag]]
      if (is.null(xcms.xcms)) next
      #xcms.se <- get_xcms_feature_se(xcms.xcms)

    }

    ### find +C +N from peaks in unlabeled sample
    {

      cn.list <- PAVE_find_xcms_CN(xcms.xcms,BPPARAM = BPPARAM)
      object@advancedAna$PAVE[[polarity.index[as.character(i.pol)]]] <- cn.list

    }



    ### timer
    {
      time.end <- Sys.time()
      time.cost <- difftime(time.end,time.start,units = "secs")%>%as.numeric()
      readr::write_lines(
        paste0("pave1",",",
               nrow(xcms::featureDefinitions(xcms.xcms)),",",
               time.cost),
        append = T,
        file = get_dir_expand_from_onedrive("Documents/YLF_Lab/Project/2025.10.10.PAVE/data/pave.CN.count.timer.csv")
      )
    }

  }

  return(object)


}


PAVE_junk_remover <- function(object,ppm = 10,rt.tol = 20){



  ### inside polarity
  {
    cn.seed.pol <- list()
    for (i.pol in 0:1) {

      pol <- ifelse(i.pol==0,"Negative","Positive")

      cn.list <- object@advancedAna$PAVE[[pol]]
      cn.seed <- lapply(cn.list,function(x){
        x %>%dplyr::mutate(pave_formula = paste0("C",max(C_count),
                                                 "N",max(N_count)))%>%
          dplyr::filter(feature_id == pave_seed)
      })%>%data.table::rbindlist()%>%
        as.data.frame()

      cn.seed <- cn.seed%>%
        dplyr::filter(pave_cor > 0.75)%>%
        dplyr::mutate(rtg = cluster_rt(rt = rtmed,rt.tol = 20),
                      pave_junkremover = "")


      ### adduct and isotope
      {

        message_with_time("Find isotope and adduct in ",pol)
        iso.diff <- data.frame(
          chemform = c("[13]C1C-1","[13]C1C-1","[18]O1O-1","[18]O1O-1","[15]N1N-1","[34]S1S-1","[37]ClCl-1"),
          charge = c(1,2,1,2,1,1,1)
        )%>%
          dplyr::mutate(
            type = "isotope",
            mass_diff = MSCC::chemform_mz(chemform)/charge )%>%
          dplyr::filter(charge%in% c(1,2))

        data("pave_adduct")
        adducts.diff <- pave_adduct%>%
          dplyr::filter(polarity == pol)%>%
          dplyr::mutate(type = "adduct")

        mass_diff <- bind_rows(iso.diff,  adducts.diff)




        cn.seed.adduct.isotope <- lapply(unique(cn.seed$rtg),
                                         function(x){

                                           #message(x)
                                           this.seed <- cn.seed %>%
                                             dplyr::filter(rtg == x)

                                           if (nrow(this.seed) <2) {
                                             return(this.seed)
                                           }

                                           this.seed.diff <-
                                             expand.grid(ft1 = seq_along(this.seed$feature_id),
                                                         ft2 = seq_along(this.seed$feature_id))%>%
                                             dplyr::filter(ft2>ft1)%>%
                                             dplyr::mutate(mz1 = this.seed$mzmed[ft1],
                                                           mz2 = this.seed$mzmed[ft2],
                                                           mzd = abs(mz1-mz2),
                                                           idx = match_mz( mzd, mass_diff$mass_diff , mz.ppm = 10000 ),
                                                           mass_diff[idx,],
                                                           mzd.ppm = abs(mzd-mass_diff)/mean(mz1+mz2)*1e6
                                             )%>%
                                             dplyr::filter(!is.na(idx),mzd.ppm < ppm)

                                           if (nrow(this.seed.diff)) {

                                             this.seed$pave_junkremover[this.seed.diff$ft2] <-
                                               this.seed.diff$type

                                           }

                                           return(this.seed)

                                         })%>%data.table::rbindlist()%>%
          as.data.frame()



        }

      ### dimer
      {

        message_with_time("Find dimer in ",pol)
        cn.seed.dimer <- lapply(unique(cn.seed$rtg),
                                function(x){

                                  #message(x)
                                  this.seed <- cn.seed %>%
                                    dplyr::filter(rtg == x)

                                  if (nrow(this.seed) <2) {
                                    return(this.seed)
                                  }

                                  this.seed.dimer <-
                                    this.seed%>%
                                    dplyr::mutate(dimer.mz = (mzmed + (i.pol-0.5)*2*1.00727) / 2 )

                                  this.seed.is.dimer <- this.seed%>%
                                    dplyr::mutate(
                                      dimer.matched = match_mz(mzmed , this.seed.dimer$dimer.mz,mz.ppm = ppm ),
                                      pave_junkremover = case_when(!is.na(dimer.matched)~ "dimer",
                                                                   T~""))

                                  this.seed$pave_junkremover <-this.seed.is.dimer$pave_junkremover


                                  return(this.seed)

                                })%>%data.table::rbindlist()%>%
          as.data.frame()



      }

      ### ring
      {


        message_with_time("Find ringing in ",pol)
        xcms.xcms <- object@xcmsData[[paste0(pol,'MS1')]]
        xcms.fdf <- xcms::featureDefinitions(xcms.xcms)%>%
          as.data.frame()
        xcms.pda <- pData(xcms.xcms)
        xcms.pave.sample <- xcms.pda%>%
          dplyr::filter(sample.type %in% c("S12C14N","S13C14N","S12C15N","S13C15N"))
        xcms.val <- xcms::featureValues(xcms.xcms, missing  = 0,value = "maxo")
        pave.sample.val <- apply(xcms.val,1,mean_f, f= xcms.pda$sample.type)%>%t
        pave.sample.val <- pave.sample.val[,c("S12C14N")]
        xcms.fdf$peakMaxo <- pave.sample.val

        cn.seed.ring <- cn.seed
        for (i in 1:nrow(cn.seed.ring)) {

          this.mz <- cn.seed.ring$mzmed[i]
          this.rt <- cn.seed.ring$rtmed[i]
          this.peakmaxo <- cn.seed.ring$peakMaxo[i]
          this.mz.range50 <- mz.range.ppm(this.mz,50)
          this.mz.range500 <- mz.range.ppm(this.mz,500)
          this.ring.maxo <- xcms.fdf %>%
            dplyr::filter(between.range(rtmed , this.rt+ c(-rt.tol,rt.tol)),
                          between.range(mzmed,c(this.mz.range500[1],this.mz.range50[1])) |
                            between.range(mzmed,c(this.mz.range50[2],this.mz.range500[2])) )%>%
            dplyr::pull(peakMaxo)%>%
            max()

          if(this.ring.maxo>this.peakmaxo*100)
            cn.seed.ring$pave_junkremover[i] <- "ringing"

        }

      }


      ### integrate
      {

        table(cn.seed.adduct.isotope$pave_junkremover)
        table(cn.seed.dimer$pave_junkremover)
        table(cn.seed.ring$pave_junkremover)

        rownames(cn.seed.adduct.isotope) <- cn.seed.adduct.isotope$feature_id
        rownames(cn.seed.dimer) <- cn.seed.dimer$feature_id
        rownames(cn.seed.ring) <- cn.seed.ring$feature_id

        cn.seed.adduct.isotope <- cn.seed.adduct.isotope[cn.seed$feature_id,]
        cn.seed.dimer <- cn.seed.dimer[cn.seed$feature_id,]
        cn.seed.ring <- cn.seed.ring[cn.seed$feature_id,]


        cn.seed$pave_junkremover <-
          integrate_anntation(cn.seed$pave_junkremover,
                              cn.seed.adduct.isotope$pave_junkremover)

        cn.seed$pave_junkremover <-
          integrate_anntation(cn.seed$pave_junkremover,
                              cn.seed.dimer$pave_junkremover)

        cn.seed$pave_junkremover <-
          integrate_anntation(cn.seed$pave_junkremover,
                              cn.seed.ring$pave_junkremover)
        cn.seed.pol[[pol]] <- cn.seed


      }




    }
  }

  ### between polarity
  {

    diff_to_neg <- pave_adduct %>%
      dplyr::filter(polarity=="Positive")%>%
      dplyr::mutate(mass_diff = mass_diff + 1.0078250320*2 - 2* 0.00054857990943 )

    diff_to_pos <- pave_adduct %>%
      dplyr::filter(polarity=="Negative")%>%
      dplyr::mutate(mass_diff = mass_diff - 1.0078250320*2 + 2* 0.00054857990943 )


    cn.seed <- data.table::rbindlist(cn.seed.pol)
    cn.seed <- cn.seed%>%
      dplyr::mutate(rtg = cluster_rt(rt = rtmed,rt.tol = 20))

    cn.seed.list <- list()
    for (i in unique(cn.seed$rtg)) {


      this.cn.seed <- cn.seed%>%
        dplyr::filter(rtg == i)

      possible.adduct.neg <- this.cn.seed %>%
        dplyr::filter(#pave_junkremover=="",
          polarity == 1)%>%
        dplyr::pull(mzmed)%>%
        expand.grid(mz = ., mzd = diff_to_pos$mass_diff)%>%
        dplyr::mutate(mz.expected = mz + mzd)


      possible.adduct.pos <- this.cn.seed %>%
        dplyr::filter(#pave_junkremover=="",
          polarity == 0)%>%
        dplyr::pull(mzmed)%>%
        expand.grid(mz = ., mzd = diff_to_neg$mass_diff)%>%
        dplyr::mutate(mz.expected = mz + mzd)

      this.cn.seed <- this.cn.seed%>%
        dplyr::ungroup()%>%
        dplyr::mutate(
          adduct.match = case_when(
            polarity == 0 ~ match_mz(mzmed,possible.adduct.neg$mz.expected,mz.ppm = ppm),
            polarity == 1 ~ match_mz(mzmed,possible.adduct.pos$mz.expected,mz.ppm = ppm)
          ),
          pave_junkremover = case_when(
            is.na(adduct.match) ~ pave_junkremover,
            T ~ integrate_anntation(pave_junkremover,"opposite_adduct")
          )
        )#%>%
      #dplyr::select(-adduct.match)
      #message(sum(!is.na(this.cn.seed$adduct.match)))

      cn.seed.list[[i]] <- this.cn.seed

    }

    cn.seed <- data.table::rbindlist(cn.seed.list)%>%
      as.data.frame()

  }


  ### Low C
  {
    data("PAVE_LowC_cutoff")
    cn.seed <- cn.seed%>%
      dplyr::mutate(
        pave_lowC_cutoff = PAVE_LowC_cutoff$mass_max[
          match(get_formula_ele_count(pave_formula,"C"),PAVE_LowC_cutoff$c.count) ],
        pave_junkremover =
          case_when(
            mzmed > pave_lowC_cutoff~integrate_anntation(pave_junkremover,"LowC"),
            T~pave_junkremover) )



  }




  ### return
  {

    for (i.pol in 0:1) {

      pol <- ifelse(i.pol==0,"Negative","Positive")
      cn.list <- object@advancedAna$PAVE[[pol]]
      cn.seed.pol <- cn.seed%>%
        dplyr::filter(polarity == i.pol)%>%
        dplyr::mutate(tmp = feature_id)%>%
        tibble::column_to_rownames("tmp")

      cn.list.junkremoved <- lapply(cn.list,function(x){

        x %>%
          dplyr::mutate(pave_formula =  cn.seed.pol[pave_seed,"pave_formula"],
                        pave_junkremover = cn.seed.pol[pave_seed,"pave_junkremover"])
      })



      cn.list.junkremoved -> object@advancedAna$PAVE[[pol]]

    }

    return(object)

  }







}

PAVE_formula_assign <- function(object,ppm = 10,rt.tol = 20){


  #object@advancedAna$PAVE
  cpdb_path = "C:/Users/91879/OneDrive/Code/R/data/MSDB/CompoundDB/CompoundDB.sqlite"
  cpdb <- CompoundDb::CompDb(cpdb_path)


  cn.list.pol <- list()
  for (i in 0:1) {

    pol <- ifelse(i==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    if (is.null(xcms.xcms)) next
    message_with_time("Find MS1 candidate...")
    xcms.xcms <- xcms_get_feature_ms1_candidate(xcms.xcms,
                                                cpdb,
                                                ppm = ppm)
    xcms.fdf <- xcms::featureDefinitions(xcms.xcms)

    cn.list <- object@advancedAna$PAVE[[pol]]
    cn.seed <- lapply(cn.list,function(x){
      x %>%dplyr::mutate(pave_formula = paste0("C",max(C_count),
                                               "N",max(N_count)))%>%
        dplyr::filter(feature_id == pave_seed)
    })%>%data.table::rbindlist()%>%
      as.data.frame()%>%
      dplyr::mutate(pave_formula_matched = F)%>%
      dplyr::filter(pave_cor > 0.75)

    for (i.seed in 1:nrow(cn.seed)) {

      this.fid <- cn.seed$feature_id[i.seed]
      this.pave.formula <- cn.seed$pave_formula[i.seed]

      this.adduct.candidate <- xcms.fdf[this.fid,]$candidate.adduct%>%unlist()
      this.chemform.candidate <- xcms.fdf[this.fid,]$candidate.formula%>%unlist()

      candidate.c.count <- get_formula_ele_count(this.chemform.candidate,"C")
      candidate.n.count <- get_formula_ele_count(this.chemform.candidate,"N")
      candidate.cn.formula <- paste0("C",candidate.c.count,"N",candidate.n.count)

      idx <- match(this.pave.formula,candidate.cn.formula)
      if(!is.na(idx))
        cn.seed$pave_formula_matched[i.seed] <- T

    }


    cn.list.pol[[pol]] <- cn.seed


  }


  ### return
  {

    for (i.pol in 0:1) {

      pol <- ifelse(i.pol==0,"Negative","Positive")
      cn.list <- object@advancedAna$PAVE[[pol]]
      cn.seed.pol <- cn.list.pol[[pol]]%>%
        dplyr::filter(polarity == i.pol)%>%
        dplyr::mutate(tmp = feature_id)%>%
        tibble::column_to_rownames("tmp")

      cn.list.formula.assigned <- lapply(cn.list,function(x){

        x %>%
          dplyr::mutate(pave_formula_matched =  cn.seed.pol[pave_seed,"pave_formula_matched"])
      })



      cn.list.formula.assigned -> object@advancedAna$PAVE[[pol]]

    }

    return(object)

  }

}



PAVE_report <- function(object,file = tempfile(fileext = "pdf"),mzr = c(0,Inf)){


  cn.stat.list <- list()
  for (i.pol in 0:1) {

    pol <- ifelse(i.pol==0,"Negative","Positive")
    xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
    xcms.fdf <- xcms::featureDefinitions(xcms.xcms)%>%
      as.data.frame()%>%
      dplyr::filter(between.range(mzmed,mzr))

    cn.list <- object@advancedAna$PAVE[[pol]]
    cn.peaks <- cn.list%>%
      data.table::rbindlist()%>%
      as.data.frame()%>%
      dplyr::filter(between.range(mzmed,mzr))%>%
      dplyr::group_by(pave_seed)%>%
      dplyr::mutate(pave_cor = na.omit(pave_cor))%>%
      dplyr::ungroup()%>%
      dplyr::group_by(feature_id)%>%
      dplyr::slice_max(pave_cor)%>%
      dplyr::ungroup()%>%
      dplyr::distinct(feature_id,.keep_all = T)

    cn.peaks.high.cor <- cn.peaks%>%
      dplyr::filter(pave_cor >= 0.75)%>%
      dplyr::mutate(
        isotope = grepl("^isotope",pave_junkremover),
        adduct = grepl("^adduct",pave_junkremover),
        LowC = grepl("^LowC",pave_junkremover),
        opposite_adduct  = grepl("^opposite_adduct",pave_junkremover),
        dimer   = grepl("^dimer",pave_junkremover),
        ringing    = grepl("^ringing",pave_junkremover)
      )

    cn.peaks.high.cor.formula <-cn.peaks.high.cor %>%
      dplyr::filter(pave_junkremover == "")




    cn.stat.list[[pol]]$ATOMCOUNT <-
      list(total_peaks = nrow(xcms.fdf),
           peaks_in_blak = 0,
           peaks_withou_labeling = length(setdiff(xcms.fdf$feature_id,cn.peaks$feature_id)),
           peaks_low_cor = sum(cn.peaks$pave_cor < 0.75),
           peaks_high_cor = sum(cn.peaks$pave_cor >= 0.75) )


    cn.stat.list[[pol]]$JUNKREMOVER <-
      list( isotopes = sum(cn.peaks.high.cor$isotope),
            adduct = sum(cn.peaks.high.cor$adduct),
            LowC = sum(cn.peaks.high.cor$LowC),
            opposite_adduct = sum(cn.peaks.high.cor$opposite_adduct),
            dimer = sum(cn.peaks.high.cor$dimer),
            ringing = sum(cn.peaks.high.cor$ringing)
      )

    cn.stat.list[[pol]]$`Formula assignment` <-
      list( formula.matched = sum(cn.peaks.high.cor.formula$pave_formula_matched),
            formula.non.matched =sum(!cn.peaks.high.cor.formula$pave_formula_matched)    )


    cn.stat.list[[pol]] <- lapply(cn.stat.list[[pol]],function(x){
      do.call(rbind,x)%>%
        `colnames<-`(pol)%>%
        as.data.frame()%>%
        tibble::rownames_to_column("PAVE annotation")
    })%>%
      rbindlist(idcol = "PAVE FUN")
  }


  ###
  {
    cn.stat.df <- cn.stat.list$Positive%>%
      dplyr::mutate(Negative = cn.stat.list$Negative$Negative)%>%
      dplyr::mutate(tmp.ratio = Positive/Positive[1],
                    tmp.ratio = num2percent(tmp.ratio),
                    #Positive = paste0(Positive,"(",tmp.ratio,")"),

                    tmp.ratio = Negative/Negative[1],
                    tmp.ratio = num2percent(tmp.ratio),
                    #Negative = paste0(Negative,"(",tmp.ratio,")"),
                    `PAVE annotation` = c(
                      "total peaks number",
                      "peaks in procedure blank",
                      "other peaks without labeling",
                      "labeling but rho<0.75",
                      "logical labeling (i.e., biological)",
                      "isotopes",
                      "dimer or double charge",
                      "adducts(assigned using same polarity mode)",
                      "adducts(assigned only using opposite polarity mode)",
                      "too low C count for mass",
                      "ringing peaks",
                      "formula match to metabolite",
                      "no formula match in database"
                    )
      )%>%
      dplyr::select(-tmp.ratio)

    # cn.stat.df%>%
    #   gt::gt(rowname_col = "PAVE annotation", groupname_col = "PAVE FUN")

    return(cn.stat.df)

  }

  ### SNR
  if(F) {

    xcms.xcms <- object@xcmsData$PositiveMS1
    xcms.xcms@.processHistory[[1]]@type <- "Peak detection"
    xcms.xcms@.processHistory[[1]]@param <- object@processingInfo$MSdevParam$findChromPeaks
    p1 <- plot_xcms_peaks_SN_distribution(xcms.xcms)

    xcms.xcms <- object@xcmsData$NegativeMS1
    xcms.xcms@.processHistory[[1]]@type <- "Peak detection"
    xcms.xcms@.processHistory[[1]]@param <- object@processingInfo$MSdevParam$findChromPeaks
    p2 <- plot_xcms_peaks_SN_distribution(xcms.xcms)
    open_plot_win(p1+p2,10,5)
  }

}

PAVE_find_xcms_CN <- function(xcms.xcms, rt.tol = 20, ppm= 10 ,
                              BPPARAM = SnowParam(workers = 6,progressbar = T) ){



  ### prepare data
  {

    CN_mass_diff_df <- get_CN_mass_diff_table(C_max = 99,N_max = 10)

    xcms.fdf <- xcms::featureDefinitions(xcms.xcms)%>%
      as.data.frame()%>%
      dplyr::mutate(pave_seed = "",
                    pave_CN = "",
                    pave_cor = NA)
    xcms.pda <- pData(xcms.xcms)
    xcms.pave.sample <- xcms.pda%>%
      dplyr::filter(sample.type %in% c("S12C14N","S13C14N","S12C15N","S13C15N"))
    xcms.val <- xcms::featureValues(xcms.xcms, missing  = 0,value = "maxo")
    #pave.sample.val <- apply(xcms.val,1,median_f, f= xcms.pda$sample.type)%>%t
    #pave.sample.val <- pave.sample.val[,c("S12C14N","S13C14N","S12C15N","S13C15N")]
    #xcms.fdf[,c("S12C14N","S13C14N","S12C15N","S13C15N")] <- pave.sample.val
  }


  ### find CN candidate
  {
    cn.list <- BiocParallel::bplapply(
      seq_along(xcms.fdf$feature_id),
      #1:300,
      FUN = function(i.ft,xcms.fdf,rt.tol,ppm,xcms.val,xcms.pave.sample){

        this.fid <- xcms.fdf$feature_id[i.ft]
        #message_with_time(this.fid)
        this.mz <- xcms.fdf$mzmed[i.ft]
        this.rt <- xcms.fdf$rtmed[i.ft]
        this.CN_mass_diff_df <- CN_mass_diff_df%>%
          dplyr::filter( C_count <= this.mz/14)

        mz.pred <- this.mz+this.CN_mass_diff_df$mass_diff
        this.fdf <- xcms.fdf%>%
          dplyr::filter(abs(rtmed - this.rt) < rt.tol,
                        mzmed >= this.mz
          )%>%
          dplyr::mutate(idx = match_mz(mzmed,mz.pred,mz.ppm = ppm),
                        this.CN_mass_diff_df[idx,]  )%>%
          dplyr::filter(!is.na(idx)   )%>%
          dplyr::mutate(
            pave_seed = paste0(pave_seed,this.fid,""),
            pave_CN = paste0("C",C_count,"N",N_count,"")
          )
        #this.fdf[,c("S12C14N","S13C14N","S12C15N","S13C15N")] <-
        #  this.fdf[,c("S12C14N","S13C14N","S12C15N","S13C15N")]/
        #  this.fdf[1,c("S12C14N")]


        possible.c.count <- unique(this.fdf$C_count)%>%setdiff(0)
        possible.n.count <- unique(this.fdf$N_count)
        cn.comb <- expand.grid(C = possible.c.count,
                               N = possible.n.count,
                               p.cor = NA)

        if (!nrow(cn.comb)) return(NULL)
        ### score possible C and N pattern
        cn.comb.list <- list()
        for (i.cn in 1:nrow(cn.comb)) {
          this.c <- cn.comb$C[i.cn]
          this.n <- cn.comb$N[i.cn]
          all.form <- c("C0N0",paste0("C0N",this.n,""),paste0("C",this.c,"N0"),paste0("C",this.c,"N",this.n,""))
          if (all(all.form %in% this.fdf$pave_CN) ) {

            cn.ft <-this.fdf%>%
              dplyr::filter(pave_CN %in% all.form)%>%
              dplyr::group_by(pave_CN)%>%
              dplyr::slice_max(peakMaxo)%>%
              dplyr::ungroup()
            cn.comb.list[[i.cn]] <- cn.ft
            m.detected <- xcms.val[cn.ft$feature_id,xcms.pave.sample$sampleNames]
            colnames(m.detected) <- xcms.pave.sample$sample.type
            rownames(m.detected) <- cn.ft$pave_CN
            m.detected <- m.detected/m.detected[1,1]
            m.ideal <- get_ideal_CN_ratio(this.c,this.n)%>%t
            m.ideal <- m.ideal[rownames(m.detected),colnames(m.detected)]

            p.cor <- cor(as.vector(m.detected),as.vector(m.ideal))
            cn.comb$p.cor[i.cn] <- p.cor

          }
        }


        ### filter xcms.fdf to return
        {
          if (any(cn.comb$p.cor > 0,na.rm =T)) {
            cn.fdf <- cn.comb.list[[which.max(cn.comb$p.cor)]]
            cn.fdf$pave_cor <- max(cn.comb$p.cor,na.rm  =T)

            #message_with_time(this.fid," Pattern: ",cn.fdf$pave_CN[-1],"; Cor = ",cn.fdf$pave_cor[1])
            return(cn.fdf)
          }

        }


        #


      },
      xcms.fdf = xcms.fdf,rt.tol =rt.tol,ppm = ppm,xcms.val=xcms.val,xcms.pave.sample=xcms.pave.sample,
      BPPARAM = BPPARAM
    )

    cn.list <- cn.list[!sapply(cn.list,is.null)]
    names(cn.list) <- sapply(cn.list,function(x) x$feature_id[1])
  }


  ### save to xcms fdf
  if(F){

    for (i.cnl in seq_along(cn.list)) {
      this.fdf <- cn.list[[i.cnl]]
      xcms.fdf[this.fdf$feature_id,]$pave_seed <-
        paste0(xcms.fdf[this.fdf$feature_id,]$pave_seed , this.fdf$pave_seed)
      xcms.fdf[this.fdf$feature_id,]$pave_CN <-
        paste0(xcms.fdf[this.fdf$feature_id,]$pave_CN , this.fdf$pave_CN)
      xcms.fdf[this.fdf$feature_id[1],]$pave_cor <- this.fdf$pave_cor[1]

    }
    xcms::featureDefinitions(xcms.xcms)$pave_seed <- xcms.fdf$pave_seed
    xcms::featureDefinitions(xcms.xcms)$pave_CN <- xcms.fdf$pave_CN
    xcms::featureDefinitions(xcms.xcms)$pave_cor <- xcms.fdf$pave_cor

  }



  return(cn.list)
}



get_CN_mass_diff_table <- function(C_max=100,N_max=20){



  ### mass and max count define
  {
    C13_mass_diff= MSCC::chemform_mz("[13]CC-1")
    N15_mass_diff= MSCC::chemform_mz("[15]NN-1")

  }

  ### mass diff matrix
  if(F){
    C_mass_diff_matrix <-
      matrix(
        rep(C13_mass_diff * (0:C_max), N_max+1),
        nrow = C_max+1
      )
    N_mass_diff_matrix <-
      matrix(
        rep(N15_mass_diff * (0:N_max), C_max+1),
        ncol = N_max+1,byrow = T
      )
    CN_mass_diff_matrix <-
      C_mass_diff_matrix+N_mass_diff_matrix
    rownames(CN_mass_diff_matrix) <- paste0("C",num2str(0:C_max))
    colnames(CN_mass_diff_matrix) <- paste0("N",num2str(0:N_max))

  }

  ### mass diff data.frame
  {

    CN_mass_diff_df <-
      expand.grid(
        C_count = 0:C_max,
        N_count = 0:N_max
      )%>%
      dplyr::mutate(
        chemform_diff = paste0("[13]C",C_count,"C-",C_count,
                               "[15]N",N_count,"N-",N_count),
        #mz = MSCC::chemform_mz(chemform_diff),
        mass_diff = C_count * C13_mass_diff + N_count * N15_mass_diff
      )
    CN_mass_diff_df <- data.table::as.data.table(CN_mass_diff_df)

  }

  return(CN_mass_diff_df)



}

get_ideal_CN_ratio <- function(C = 10 , N = 2, ratio.adjust = c(1, 1, 1, 1)){

  m <- diag(rep(1,4))
  colnames(m) <- c("C0N0",paste0("C0N",N),paste0("C",C,"N0"),paste0("C",C,"N",N))
  rownames(m) <- c("S12C14N","S12C15N","S13C14N","S13C15N")

  if (length(ratio.adjust) != 4) {
    stop("ratio.adjust must have length 4.")
  }
  m <- m * as.numeric(ratio.adjust)

  if(N==0) m <- m[,c(1,3)]+m[,c(2,4)]
  m <- rbind(m,Blank = 0)
  return(m)
}


get_PAVE_LowC_cutoff <- function( c_max = 100){

  hmdb.cp <- MSdb:::get_HMDB_Compound_DF()
  atom.count <- MSCC:::chemform_parse(hmdb.cp$chemform)

  hmdb.pave.stat <- data.frame(
    chemform = hmdb.cp$chemical_formula,
    mass = as.numeric(hmdb.cp$monisotopic_molecular_weight),
    c.count =atom.count[,"C"]
  )%>%
    dplyr::filter(!is.na(mass),!is.na(c.count))


  PAVE_LowC_cutoff <- data.frame(
    c.count = 1:c_max,
    mass_min = NA,
    mass_max = NA
  )
  for (i.c in 1:c_max) {

    c.count.cp <- hmdb.pave.stat%>%
      dplyr::filter(c.count >= i.c-2,c.count <= i.c +2)

    c.c.m.q <- quantile(c.count.cp$mass,c(0.01,0.99))
    PAVE_LowC_cutoff$mass_min[i.c] <- c.c.m.q[1]
    PAVE_LowC_cutoff$mass_max[i.c] <- c.c.m.q[2]
  }

  return(PAVE_LowC_cutoff)

}

integrate_anntation <- function(str1,str2,sep = ";"){

  str_sep <- ifelse(str1==""|is.na(str1)|str2 == "","",sep)

  paste0(str1,str_sep,str2)

}


get_pave_ig_vertex_form <- function(ring.ig){


  v.form <- lapply(names(V(ring.ig)),function(x){

    #message(x)
    e <- igraph::incident(ring.ig,x,mode = "all")
    e.ends <- igraph::ends(ring.ig,e)
    ef <- e[e.ends[,1] == x]
    et <- e[e.ends[,2] == x]

    e.iso <- et$element%>%na.omit()%>%unique()
    if (isEmpty(e.iso)) e.iso <- "Natural"
    e.adduct <- c(ef$adduct.from,et$adduct.to)%>%na.omit()%>%unique()
    e.frag <- ef$fragment%>%na.omit()
    if (length(e.frag) > 0) e.frag <- MSCC::chemform_sum(e.frag)
    paste0(e.iso,";",e.adduct,";",e.frag)

  })
  names(v.form) <- names(V(ring.ig))
  return(v.form)
}

ring.node.form.group <- function(x){

  if (nrow(x)==1) return(setNames(1,1))
  xcomb <- combn(1:nrow(x),2)
  xcomb.eq <- rep(F,ncol(xcomb))
  for (i in 1:ncol(xcomb)) {
    f1 <- x[xcomb[1,i],]%>%str_split(pattern = ";")%>%lapply(function(y){
      y[y==""]<-NA
      y})
    f2 <- x[xcomb[2,i],]%>%str_split(pattern = ";")%>%lapply(function(y){
      y[y==""]<-NA
      y})
    xcomb.eq[i] <- sapply(seq_along(f1),function(z){
      all(f1[[z]] == f2[[z]],na.rm = T)
    })%>%all(na.rm = T)
    #xcomb.eq[i] <-all(x[xcomb[1,i],]==x[xcomb[2,i],],na.rm = T)

  }
  x.group <- get_igraph_membership(igraph::graph_from_data_frame(
    t(xcomb[,xcomb.eq]),vertices = 1:nrow(x)))
  return(x.group)
}






chemform_simplify <- function(chemform){

  ele.matrix <- MSCC::chemform_parse(chemform)
  ele.matrix <- ele.matrix[,order(colnames(ele.matrix))]
  MSCC:::chemform_from_ele_matrix(ele.matrix)

}


chemform_remove_iso <- function(chemform){

  ele.matrix <- MSCC::chemform_parse(chemform)
  iso.idx <- colnames(ele.matrix)[MSdev:::is.isotope(colnames(ele.matrix))]
  if (length(iso.idx)) {
    base.idx <- MSdev:::get_ele_uniso(iso.idx)
    # chemform_parse() returns a matrix, which cannot create new columns by name,
    # so add any missing base-element columns (e.g. a pure "[13]C" diff has no "C").
    missing.base <- setdiff(unique(base.idx), colnames(ele.matrix))
    if (length(missing.base)) {
      add <- matrix(
        0, nrow = nrow(ele.matrix), ncol = length(missing.base),
        dimnames = list(rownames(ele.matrix), missing.base)
      )
      ele.matrix <- cbind(ele.matrix, add)
    }
    # accumulate per isotope so duplicate base elements (e.g. [13]C and [14]C) sum
    for (k in seq_along(iso.idx)) {
      ele.matrix[, base.idx[k]] <- ele.matrix[, base.idx[k]] + ele.matrix[, iso.idx[k]]
    }
    ele.matrix <- ele.matrix[, !MSdev:::is.isotope(colnames(ele.matrix)), drop = FALSE]
  }
  MSCC:::chemform_from_ele_matrix(ele.matrix)

}





get_adduct_mass_diff <- function(polarity = 0,direction = 1){


  pol <- ifelse(polarity==1,"positive","negative")

  adduct.table <- MSCC::adduct.table |>
    dplyr::filter(
      Ion_mode == pol,
      Multi == 1,
      abs(Charge) == 1
    ) |>
    dplyr::mutate(m_c = Multi / Charge)


  adduct.diff <- expand.grid(
    adduct.from = 1:nrow(adduct.table),
    adduct.to = 1:nrow(adduct.table)
  ) |>
    dplyr::filter(
      adduct.table$m_c[adduct.from] == adduct.table$m_c[adduct.to]
    ) |>
    dplyr::mutate(
      chemform_diff = MSCC::chemform_calc(adduct.table$Formula_diff[adduct.to],
                                          adduct.table$Formula_diff[adduct.from],
                                          calc = "-",return = "chemform"),
      chemform_diff = chemform_simplify(chemform_diff),
      mass_diff = MSCC::chemform_mz(chemform_diff),
      #mass_diff =adduct.table$Mass[adduct.to] - adduct.table$Mass[adduct.from],
      #charge = adduct.table$Charge[adduct.to],
      adduct.from = adduct.table$Adduct[adduct.from],
      adduct.to = adduct.table$Adduct[adduct.to]
    )

  #which(upper.tri(diag(nrow(adduct.table)),diag = F),arr.ind = T)
  adduct.diff <- data.table::as.data.table(adduct.diff)

  return(adduct.diff)

}


get_iso_mass_diff <- function(){


  iso.ele <- c("[13]C", "[2]H", "[18]O","[15]N","[34]S",
               #"[41]K","[44]Ca","[10]B","[29]Si","[30]Si","[53]Cr", "[60]Ni","[62]Ni"
               "[37]Cl","[81]Br"
  )

  ele <- MSCC::elem_table |>
    dplyr::mutate(ele.base = MSdev:::get_ele_uniso(element)) |>
    dplyr::group_by(ele.base) |>
    dplyr::filter(any(element %in% iso.ele)) |>
    dplyr::arrange(mass) |>
    dplyr::mutate(
      chemform_diff = paste0(element, "1", element[1], "-1"),
      mass_diff = MSCC::chemform_mz(chemform_diff)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(mass_diff != 0) |>
    dplyr::select("element", "chemform_diff", "mass_diff")

  data.table::as.data.table(ele)


}


#' @import data.table
get_fragment_mass_diff <- function(){

  fragment.formula <- c(
    "C1O2","C1H2O1","H2O1","N1H3",### from netID
    "H1N1O1","C1H2","C2H4","C3H6", ### from data,
    ### 1. Basic small molecules and polar groups ###
    "H2O1",      # Water loss
    "N1H3",      # Ammonia loss
    "C1O2",      # Carbon dioxide loss (Decarboxylation)
    "C1H2O2",    # Formic acid loss
    "C1O1",      # Carbon monoxide loss

    ### 2. Lipid and energy metabolism headgroups ###
    "C2H8N1O4P1",# Phosphoethanolamine loss (PE lipid headgroup)
    "C3H5N1O2",  # Serine residue loss (PS lipid headgroup)
    "H3O4P1",    # Phosphoric acid loss (often complementary to H2PO4- in MS)

    ### 3. Sugars and Phase II metabolism modifications ###
    "C6H10O5",   # Hexose loss (e.g., glucose, galactose)
    "C5H8O4",    # Pentose loss (e.g., ribose, xylose)
    "C6H10O4",   # Deoxyhexose loss (e.g., rhamnose, fucose)
    "C6H8O6",    # Glucuronic acid loss (hallmark of glucuronidation)
    "S1O3",      # Sulfur trioxide loss (Desulfation)
    "C12H20O10", # Disaccharide loss (e.g., rutinose in flavonoids)

    ### 4. Common modifications: Methylation & Acetylation ###
    "C1H4O1",    # Methanol loss (common in methyl esters or methoxy cleavage)
    "C2H4O2",    # Acetic acid loss (classic marker for acetylated metabolites)
    "C1H2O1",    # Formaldehyde loss (common in methoxy-containing aromatics)
    "C2H2O1",    # Ketene loss (alternative form of acetyl skeleton cleavage)
    "C3H2O3",    # Malonyl group loss (common acylation in plant metabolites)

    ### 5. Nitrogen/Sulfur-specific losses ###
    "H2S1",      # Hydrogen sulfide loss (indicates sulfur-containing amino acids)
    "C1H1N1",    # Hydrogen cyanide loss (alkaloids or N-containing heterocycles)
    "C1H1N1O1",  # Isocyanic acid loss (urea metabolites or pyrimidine cleavage)
    "C1H3N1O1",  # Formamide loss (deep cleavage of N-containing rings)
    "C1H5N1",    # Methylamine loss (N-methylated metabolites)

    ### 6. Alkyl chain and hydrocarbon skeleton cleavages ###
    "C2H4",      # Ethylene loss (McLafferty rearrangement or long-chain lipid cleavage)
    "C3H6",      # Propylene loss (isoprenoid or lipid cleavage)
    "C4H8"       # Butene loss

  )

  data.table::data.table(chemform_diff = unique(fragment.formula))[,mass_diff := MSCC::chemform_mz(chemform_diff)][
    ,chemform_diff := chemform_simplify(chemform_diff)
  ][
    ,fragment := chemform_diff ]

}


pave_igraph_contract <- function(pave.ig){

  eda <- edata(pave.ig) %>%
    dplyr::filter(!type == "adduct")
  vm <- get_igraph_membership(graph_from_data_frame(eda))
  vda <- vdata(pave.ig)

  vm2 <- setNames(seq_along(vda$name)+length(vda$name),vda$name)
  vm2[names(vm)] <-vm

  ign <-igraph::contract(pave.ig,
                         as.numeric(factor(vm2)),
                         vertex.attr.comb = "first")
  #vis_pave_igraph(ign)
  return(ign)
}

