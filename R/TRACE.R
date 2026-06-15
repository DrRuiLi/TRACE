TRACE <- function(object){


  if(F){
    object <- readRDS("temp/20251119.rds")
    #object <- MSdev_load("d:/data/2026.01.07.TRACE.Nutrition/MSdev_2026_01_07.Rdata")

    demo.file <- get_dir_expand_from_onedrive(
      "Documents/YLF_Lab/Project/2025.10.10.TRACE/data/demo/TRACE.demo.rdata")
    demo.file <- "d:/data/2025.12.26.PAVE2/PAVE_With_Params/QE_Plus_ppm10_sn10.rdata"
    object <- MSdev_load(demo.file)
    object <- MSdev_set_param(
      object,
      findChromPeaks =
        xcms::CentWavePredIsoParam(
          ppm = ppm,
          prefilter = c(3,1000),
          peakwidth = c(10,30),
          snthresh = 10,
          fitgauss = T),
      groupChromPeaks =
        xcms::PeakDensityParam(
          sampleGroups = "A",
          minFraction = 0.6,
          binSize = 0.002,
          bw = 5,
          ppm = 5))
    object <- MSdev_xcmsProcessing(object)
    i.pol = 0
  }


  for (i.pol in 0:1) {

    rt.tol = 10
    ppm = 5

    ### data
    {

      time.start <- Sys.time()
      pol <- ifelse(i.pol==0,"Negative","Positive")
      message_with_time(pol)

      xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
      #xcms.xcms <- xcms_filter_feature_mz_rsd(xcms.xcms,rsd.ppm = 1.2)
      #xcms.xcms <- xcms_filter_feature_rt_rsd(xcms.xcms,rt.shift = 3.5)
      xcms.xcms <- xcms_get_feature_wmean(xcms.xcms)

      xcms.net <- MSdev:::get_xcms_feature_connect(xcms.xcms,rt.tol = rt.tol)
      xcms.fdf <- xcms::featureDefinitions(xcms.xcms)
      xcms.val <- xcms::featureValues(xcms.xcms,missing = 0,value = "maxo")
      xcms.TRACE.sample <- Biobase::pData(xcms.xcms)%>%
        dplyr::filter(sample.type %in% c("S12C14N","S12C15N","S13C14N","S13C15N","Blank"))

    }

    ### theoritical mass diff match
    {

      cn.mass.diff <- get_CN_mass_diff_table(N_max = 10)[,type :="CN_label"]
      ad.mass.diff <- get_adduct_mass_diff(polarity = i.pol)[,type := "adduct"]
      is.mass.diff <- get_iso_mass_diff()[,type := "isotope"]
      fg.mass.diff <- get_fragment_mass_diff()[,type := "fragment"]

      mass.diff.range <- range(cn.mass.diff$mass_diff,
                               ad.mass.diff$mass_diff,
                               is.mass.diff$mass_diff)
      xcms.net <- xcms.net[between(mz.diff,  mass.diff.range[1],mass.diff.range[2])]

      cn.match <- MSdev:::match_mz_foverlaps( xcms.net$mz.diff,cn.mass.diff$mass_diff,
                                      ppm.base = xcms.net$mz.mean,ppm = ppm)

      ad.match <- MSdev:::match_mz_foverlaps( xcms.net$mz.diff,ad.mass.diff$mass_diff,
                                      ppm.base = xcms.net$mz.mean,ppm = ppm)

      is.match <- MSdev:::match_mz_foverlaps( xcms.net$mz.diff,is.mass.diff$mass_diff,
                                      ppm.base = xcms.net$mz.mean,ppm = ppm)

      fg.match <- MSdev:::match_mz_foverlaps( xcms.net$mz.diff,fg.mass.diff$mass_diff,
                                      ppm.base = xcms.net$mz.mean,ppm = ppm)


      cn.net <- cbind(xcms.net[cn.match$ion1,],cn.match[,c("mz.ppm","ion1") ],
                      cn.mass.diff[cn.match$ion2,])[mass_diff > 0]
      ad.net <- cbind(xcms.net[ad.match$ion1,],ad.match[,c("mz.ppm","ion1") ],
                      ad.mass.diff[ad.match$ion2,])[mass_diff > 0]
      is.net <- cbind(xcms.net[is.match$ion1,],is.match[,c("mz.ppm","ion1") ],
                      is.mass.diff[is.match$ion2,])[mass_diff > 0]
      fg.net <- cbind(xcms.net[fg.match$ion1,],fg.match[,c("mz.ppm","ion1") ],
                      fg.mass.diff[fg.match$ion2,])[mass_diff > 0]

    }

    ### filter CN label pattern
    {

      cn.net <- cn.net%>%
        dplyr::mutate(TRACE_pattern = paste0("C",C_count ,"N",N_count  ))
      cn.net.list <- split(cn.net,cn.net$from)
      prefilt <- sapply(cn.net.list,function(x.cn){ 0 %in% x.cn$N_count   })
      cn.net.list <- cn.net.list[prefilt]
      message_with_time("Find CN label pattern...")
      cn.net.list.hit <- bplapply(names(cn.net.list),function(x){

        #message(x)
        x.cn <- cn.net.list[[x]]
        possible.c.count <- unique(x.cn$C_count)
        possible.n.count <- unique(x.cn$N_count)
        c.max <- x.cn$from.mz[1]/14
        #possible.c.count <- possible.c.count[possible.c.count < c.max&possible.c.count > 0]

        possible.c.count <- TRACE_LowC_cutoff%>%
          dplyr::filter(mass_min < x.cn$from.mz[1],
                        mass_max > x.cn$from.mz[1])%>%
          dplyr::pull(c.count)%>%intersect(possible.c.count,.)

        possible.n.count <- possible.n.count[possible.n.count < c.max]
        cn.comb <- expand.grid(C = possible.c.count,
                               N = possible.n.count,
                               p.cor = NA)

        if (!nrow(cn.comb)) return(NULL)
        for (i.cn in 1:nrow(cn.comb)) {

          this.c <- cn.comb$C[i.cn]
          this.n <- cn.comb$N[i.cn]
          all.form <- c(paste0("C0N",this.n,""),paste0("C",this.c,"N0"),
                        paste0("C",this.c,"N",this.n,""))
          all.form <- setdiff(all.form,"C0N0")
          if (!all(all.form %in% x.cn$TRACE_pattern ) ) next

          #message(x)
          to.id <- x.cn$to[match(all.form,x.cn$TRACE_pattern)]
          m.detected <- xcms.val[c(x.cn$from[1],to.id),  xcms.TRACE.sample$sampleNames]
          colnames(m.detected) <- xcms.TRACE.sample$sample.type
          rownames(m.detected) <- c("C0N0",all.form)
          mean.c0n0 <- mean(m.detected[rownames(m.detected) == "C0N0",
                                       colnames(m.detected) == "S12C14N"])
          m.detected <- m.detected/mean.c0n0
          m.ideal <- get_ideal_CN_ratio(this.c,this.n)%>%t
          m.ideal <- m.ideal[rownames(m.detected),colnames(m.detected)]
          p.cor <- cor(as.vector(m.detected),as.vector(m.ideal))
          #p.cor <- weights::wtd.cor(as.vector(m.detected),as.vector(m.ideal),weight = as.vector(m.ideal)+0.1)[1]
          cn.comb$p.cor[i.cn] <- p.cor
          if(F){
            ComplexHeatmap::Heatmap(
              cbind(m.detected,m.ideal),
              col = colramp(),
              cluster_rows = F,
              cluster_columns = F
            )
          }
        }

        p.cor.max <- max(cn.comb$p.cor,na.rm = T)
        #message(p.cor.max)
        if(p.cor.max< 0) return(NULL)
        cn.comb <- cn.comb%>%dplyr::slice_max(p.cor,with_ties = F)
        all.form <- c(paste0("C0N",cn.comb$N,""),paste0("C",cn.comb$C,"N0"),
                      paste0("C",cn.comb$C,"N",cn.comb$N,""))
        all.form <- setdiff(all.form,"C0N0")
        x.cn <- x.cn[match(all.form,x.cn$TRACE_pattern),]
        x.cn$TRACE_cor <- p.cor.max
        x.cn$TRACE_formula <-  paste0("C",cn.comb$C,"N",cn.comb$N,"")
        return(x.cn)

      },BPPARAM = SerialParam(progressbar = T))
      names(cn.net.list.hit) <- names(cn.net.list)
      cn.net.list.hit <- cn.net.list.hit[!sapply(cn.net.list.hit,is.null)]
      cn.net.hit <- data.table::rbindlist(cn.net.list.hit)%>%
        dplyr::filter(TRACE_cor > 0.75)




      ### timer and temp save
      {
        time.end <- Sys.time()
        time.cost <- difftime(time.end,time.start,units = "secs")%>%as.numeric()
        readr::write_lines(
          paste0("TRACE2",",",
                 nrow(featureDefinitions(xcms.xcms)),",",
                 time.cost),
          append = T,
          file = get_dir_expand_from_onedrive("Documents/YLF_Lab/Project/2025.10.10.TRACE/data/TRACE.CN.count.timer.csv")
        )

        if (T) {

          cn.temp <- data.table::rbindlist(cn.net.list.hit)
          object@advancedAna$TRACE2_temp[[pol]]$CNfinder <- cn.net %>%
            dplyr::mutate(TRACE_cor = cn.temp$TRACE_cor[match(ion1,cn.temp$ion1)])

        }

        }
    }


    ### RT and mz error evaluation
    if(T){



      cn.net.eval <- cn.net%>%
        dplyr::mutate(cn.hit = ion1%in% cn.net.hit$ion1)%>%
        dplyr::arrange(cn.hit)
      data.table::setDT(cn.net.eval)
      cn.net.eval[(cn.hit),mz.ppm]
      #cn.net.eval <- cn.net.eval[1:1000000,]
      ppm.fit <- distinct_norm_from_random_backgroud(cn.net.eval[(cn.hit),mz.ppm],
                                                     cn.net.eval[!(cn.hit),mz.ppm])
      ppm.dyn <- ppm.fit$sd * qnorm(0.999)
      rt.fit <- distinct_norm_from_random_backgroud(cn.net.eval[(cn.hit),rt.diff],
                                                    cn.net.eval[!(cn.hit),rt.diff])
      rt.tol.dyn <- rt.fit$sd * qnorm(0.99999)

      object@advancedAna$TRACE2_temp[[pol]][["mz.dyn"]] <- ppm.fit
      object@advancedAna$TRACE2_temp[[pol]][["rt.dyn"]] <- rt.fit

      cols <- c("TRUE" = "red","FALSE" = "#888888")

      p <- ggplot() +
        ggrastr::rasterise(
          geom_point(data = cn.net.eval,
                     aes(x = mz.ppm, y = rt.diff,
                         col = cn.hit),
                     pch = 16,alpha = 0.2,size = 0.02),
          dpi = 300)+
        scale_color_manual(values = cols)+
        labs(x = "mz error (ppm)",y = "rt shift (s)")+
        #coord_fixed(ppm/rt.tol)+
        theme_bw(base_size = 6)+
        theme(legend.position = "none")
      p.r <- ggplot(cn.net.eval)+
        geom_histogram(aes(y = rt.diff,x = after_stat(density), fill = cn.hit),
                       position = "dodge",#stat = "density",
                       bins = 20,col = "white")+
        stat_ecdf(aes(y = rt.diff, col = cn.hit),linewidth = 0.5)+
        geom_hline(yintercept = rt.tol.dyn*c(-1,1),col = "red",linewidth = 0.5 ,lty = "dashed")+
        annotate(geom = "text",x = 0.5, y = rt.tol/2,label = str_digit(rt.tol.dyn),
                 size = 2,col = "red",check_overlap = T)+
        scale_fill_manual(values =cols)+
        scale_color_manual(values =cols)+
        scale_x_continuous(expand = c(0,0),breaks = c(0,1))+
        labs(x = NULL, y = NULL)+
        theme_classic()+
        theme_classic(base_size = 6)+
        theme(axis.text.y = element_blank(),
              legend.position = "inside",
              legend.position.inside = c(0.5,0.9),
              axis.ticks = element_blank(),
              legend.title = element_text(size = 5,face = "bold"),
              legend.text = element_text(size = 4),
              legend.background = element_blank(),
              legend.key.size =unit(0.1,"inch"),
              legend.key.spacing = unit(0.02,"inch"),
              legend.title.position = "top")

      p.u <- ggplot(cn.net.eval)+
        geom_histogram(aes(x = mz.ppm,y = after_stat(density), fill = cn.hit),
                       position = "dodge",#stat = "density",
                       bins = 20,col = "white")+
        stat_ecdf(aes(x = mz.ppm, col = cn.hit),linewidth = 0.5,show.legend = F)+
        geom_vline(xintercept = ppm.dyn*c(-1,1),col = "red",linewidth = 0.5 ,lty = "dashed")+
        annotate(geom = "text",x = -ppm/2, y = 0.8,label = str_digit(ppm.dyn),
                 size = 2,col = "red",check_overlap = T)+
        scale_fill_manual(values = cols)+
        scale_color_manual(values =cols)+
        scale_y_continuous(expand = c(0,0),breaks = c(0,1))+
        labs(x = NULL, y = NULL,fill = "CN labeled")+
        theme_classic(base_size = 6)+
        theme( legend.position = "none",
               axis.text.x = element_blank(),
               axis.ticks = element_blank())

      p.ur <- ggplot(cn.net.eval)+
        geom_bar(aes( y =  0, fill = cn.hit),position = "stack")+
        annotate(geom = "text",x = 0, y = 0,size = 2,
                 label = num2percent(sum(cn.net.eval$cn.hit)/nrow(cn.net.eval)))+
        scale_fill_manual(values =cols,guide = guide_legend(ncol = 1))+
        labs(fill = "CN labeled")+
        coord_polar()+
        theme_void(base_size = 6)+
        theme(legend.title = element_text(size = 5,face = "bold"),
              legend.text = element_text(size = 5),
              legend.key.size =unit(0.1,"inch"),
              legend.key.spacing = unit(0.02,"inch"),
              legend.position = "none",
              legend.title.position = "top")
      #p.ur
      p.all <- p.u+p.ur+p+p.r+
        plot_layout(heights  = c(0.2,0.8),widths = c(0.8,0.2))+
        plot_annotation(title = paste0(get_MSdev_instrument(object)," ",pol))&
        theme(
          plot.title = element_text(size = 6),
          plot.tag.position = "topleft",
          plot.margin = margin(t = 1, r = 1 , b = 1, l = 1 ))
      #open_plot_win(p.all,width = 10,height = 10)
      fo <- paste0(object@projectInfo$MSdevFile,".TRACE.error.pdf")
      fo <- paste0( "C:\\Users\\91879\\OneDrive\\Documents\\YLF_Lab\\Project\\2025.10.10.TRACE\\result/dynamic error/",
                    basename(fo) )
      export_graph2pdf(p.all,file_path = fo,width = 3,height = 3,append = i.pol)




    }


    ### filter with dynamic error range
    {


      cn.net.filter <- cn.net.hit[abs(mz.ppm )> ppm.dyn | abs(rt.diff) > rt.tol.dyn]
      cn.net.hit <- cn.net.hit[! from %in% cn.net.filter$from]
      ad.net <- ad.net[abs(mz.ppm )< ppm.dyn & abs(rt.diff) < rt.tol.dyn]
      is.net <- is.net[abs(mz.ppm) < ppm.dyn & abs(rt.diff) < rt.tol.dyn]
      fg.net <- fg.net[abs(mz.ppm) < ppm.dyn & abs(rt.diff) < rt.tol.dyn]

    }


    ### CN seed network
    {

      cn.seed <- cn.net.hit$from
      cn.seed.formula <- cn.net.hit[,.SD[1],by = from]
      cn.seed.formula <- setNames(cn.seed.formula$TRACE_formula,cn.seed.formula$from)
      ad.net.cs <- ad.net[from %in% cn.seed&to %in% cn.seed]
      is.net.cs <- is.net[from %in% cn.seed&to %in% cn.seed]
      fg.net.cs <- fg.net[from %in% cn.seed&to %in% cn.seed]
      cn.seed.net <- bind_rows(ad.net.cs,is.net.cs,fg.net.cs) %>%
        data.table::as.data.table()

      .equal.cn.chemform <- function(cn.diff,chemfrom.diff){
        m <- MSCC::chemform_parse(c(cn.diff,chemfrom.diff))
        m <- get_matrix_value_fill_with_NA(m,colnames_vec = c("C","N"))
        m[is.na(m)] <- 0
        m1 <- m[seq_along(cn.diff),]
        m2 <- m[length(cn.diff)+seq_along(cn.diff),]
        eq <- m1[,"C"]==m2[,"C"]&m1[,"N"]==m2[,"N"]
        unname(eq)
      }
      cn.seed.net[
        ,chemform_diff := chemform_simplify(chemform_diff)][
          ,temp := data.table::fcase(
            type == "adduct",paste0(adduct.from,">>",adduct.to),
            default = chemform_diff)
        ][,label :=  paste0(type,": ",temp) ][
          ,from.cn := cn.seed.formula[as.character(from)]] [,
                                                            to.cn := cn.seed.formula[as.character(to)]][,
                                                                                                        cn.diff := MSCC::chemform_calc(to.cn,from.cn ,"-",return = "chemform")][
                                                                                                          ,chemform.equal := .equal.cn.chemform(cn.diff,chemform_diff)]
      cn.seed.net <- cn.seed.net %>%
        dplyr::mutate(
          chemform_diff := chemform_simplify(chemform_diff),
          label = case_when(
            type == "adduct"~paste0(adduct.from,">>",adduct.to),
            T~ chemform_diff
          ),
          label = paste0(type,": ",label),
          from.cn = cn.seed.formula[as.character(from)] ,
          to.cn = cn.seed.formula[as.character(to)] ,
          cn.diff = MSCC::chemform_calc(to.cn,from.cn ,"-",return = "chemform"),
          cn.equal = cn.diff == "",
          chemform.equal = .equal.cn.chemform(cn.diff,chemform_diff),
          new.type = case_when(
            type == "adduct"&cn.equal ~ "adduct",
            #type == "adduct"&!cn.equal&chemform.equal ~ "fragment",
            type == "adduct"&!cn.equal&chemform.equal ~ "false",
            type == "adduct"&!cn.equal&!chemform.equal ~ "false",

            type == "fragment"&chemform.equal ~ "fragment",
            type == "fragment"&!chemform.equal ~ "false",

            type == "isotope"&chemform.equal ~ "isotope",
            type == "isotope"&!chemform.equal&element=="[13]C"&cn.diff == "C-2" ~ "isotope",

            T~"false"
          ),
          old.type = type,
          type = new.type
        )%>%
        dplyr::filter(type != "false")%>%
        dplyr::mutate(eid = 1:n())
      cn.seed.node <- data.frame(
        name = as.character(unique(cn.seed))
      )%>%
        dplyr::mutate(
          color = case_when(
            name %in% cn.net.hit$from ~ "#E64B35",
            T~"#97C2FC"  ),
          TRACE_formula = cn.seed.formula[name],
          mz = xcms.fdf$mzmed[as.numeric(name)],
          rt = xcms.fdf$rtmed[as.numeric(name)],
          label = paste0(TRACE_formula,"\n",name))
      cn.seed.ig <- igraph::graph_from_data_frame(cn.seed.net,vertices = cn.seed.node)


    }

    ### CN seed candidate
    {

      cn.seed.vdata <- vdata(cn.seed.ig)
      cpdb <- openxlsx::read.xlsx("d:/data/2025.12.26.TRACE2/trace.cp.db.xlsx")
      adducts <- MSCC::adduct.table%>%
        dplyr::filter((sign(Charge)+1)/2 == unique(polarity(xcms.xcms)),
                      Multi  == 1,
                      abs(Charge) == 1)
      cp.adduct <- MSCC::chemform_adduct(cpdb$formula,
                                         adducts$Adduct,
                                         value = "all" )
      cp.adduct <- cp.adduct%>%
        dplyr::mutate(compound_id= cpdb$compound_id[id] ,
                      rt = cpdb$rt[id] )%>%
        dplyr::filter( findInterval(chemform.adduct.mz,
                                    mzrange(xcms.xcms))==1)
      matched.df <- match_mz_foverlaps(mz1 = cn.seed.vdata$mz,
                                       mz2 = cp.adduct$chemform.adduct.mz,
                                       ppm = 10)
      matched.df2 <- cbind( matched.df,cp.adduct[matched.df$ion2,])
      cn.seed.vdata$candidate.id <- sapply(1:nrow(cn.seed.vdata),function(i){
        idx <- matched.df$ion2[matched.df$ion1 == i]
        cp.adduct$compound_id[as.numeric(idx)]
      })
      cn.seed.vdata$candidate.formula <- sapply(1:nrow(cn.seed.vdata),function(i){
        idx <- matched.df$ion2[matched.df$ion1 == i]
        cp.adduct$chemform[as.numeric(idx)]
      })
      cn.seed.vdata$candidate.adduct <- sapply(1:nrow(cn.seed.vdata),function(i){
        idx <- matched.df$ion2[matched.df$ion1 == i]
        cp.adduct$adduct[as.numeric(idx)]
      })
      cn.seed.vdata$candidate.rt <- sapply(1:nrow(cn.seed.vdata),function(i){
        idx <- matched.df$ion2[matched.df$ion1 == i]
        cp.adduct$rt[as.numeric(idx)]
      })

      for (x in 1:nrow(cn.seed.vdata))  {

        x.name <- cn.seed.vdata$name[x]
        x.cn <- cn.seed.vdata$TRACE_formula[x]
        x.candi.id <- cn.seed.vdata$candidate.id[[x]]
        x.candi.formula <- cn.seed.vdata$candidate.formula[[x]]
        x.candi.cn <-  extract_formula_CN(x.candi.formula)
        x.candi.adduct <- cn.seed.vdata$candidate.adduct[[x]]
        x.candi.rt <- cn.seed.vdata$candidate.rt[[x]]

        id.cn.match <- x.candi.cn %in% x.cn
        x.adduct <- x.candi.adduct[id.cn.match]


        x.candi.id[id.cn.match] -> cn.seed.vdata$candidate.id[[x]]
        x.candi.formula[id.cn.match] -> cn.seed.vdata$candidate.formula[[x]]
        x.candi.adduct[id.cn.match] -> cn.seed.vdata$candidate.adduct[[x]]
        x.candi.rt[id.cn.match] -> cn.seed.vdata$candidate.rt[[x]]
      }
    }


    ### TRACE annotate
    if(T){


      node.annos <- lapply(1:nrow(cn.seed.vdata), function(x){

        x.name <- cn.seed.vdata$name[x]
        x.cn <- cn.seed.vdata$TRACE_formula[x]
        x.candi.id <- cn.seed.vdata$candidate.id[[x]]
        x.candi.formula <- cn.seed.vdata$candidate.formula[[x]]
        x.candi.adduct <- cn.seed.vdata$candidate.adduct[[x]]
        x.candi.rt <- cn.seed.vdata$candidate.rt[[x]]
        x.candi.rt[is.na(x.candi.rt)] <- Inf
        x.candi.rtd <- x.candi.rt - cn.seed.vdata$rt[x]

        x.from <- edata(cn.seed.ig)%>%
          dplyr::filter(type != "isotope",
                        from == x.name)%>%
          dplyr::select(type,eid,
                        adduct = adduct.from,
                        seed = to,
                        fragment,element
          )

        x.to <- edata(cn.seed.ig)%>%
          dplyr::filter(type != "fragment",
                        to == x.name ) %>%
          dplyr::select(type,eid,
                        adduct =adduct.to,
                        seed = from,
                        fragment,element
          )

        anno <- bind_rows(x.from,  x.to)

        ad.score <- ifelse(x.candi.adduct%in% anno$adduct, 1, 0)
        rt.score <- (1 - abs(x.candi.rtd)/1000 )
        rt.score[rt.score < 0] <- 0
        #rt.score <- sign(rt.score)
        score <- (ad.score + rt.score) / 2

        if(any((score)) > 0 ){

          idx <- which.max(score)
          res <- data.frame(
            name = x.name,
            type = "metabolite",
            seed = x.name,
            score = score[idx],
            compound.id = x.candi.id[idx],
            compound.formula = x.candi.formula[idx],
            compound.adduct = x.candi.adduct[idx],
            compound.rt = x.candi.rt[idx]
          )

        }else if( any(c("fragment","isotope") %in% anno$type ) ){

          x.anno <- anno %>%
            dplyr::filter(type %in% c("fragment","isotope"))%>%
            dplyr::slice_head(n=1)

          res <- data.frame(
            name = x.name,
            type = x.anno$type,
            seed = x.anno$seed,
            score = 0,
            compound.id = NA,
            compound.formula = NA,
            compound.adduct = NA,
            compound.rt = NA
          )


        }else if (length(score) > 0 ){

          res <- data.frame(
            name = x.name,
            type = "metabolite",
            seed = x.name,
            score = 0,
            compound.id = x.candi.id[1],
            compound.formula = x.candi.formula[1],
            compound.adduct = x.candi.adduct[1],
            compound.rt = x.candi.rt[1]
          )

        }else{

          res <- data.frame(
            name = x.name,
            type = "unknown",
            seed = NA,
            score = 0,
            compound.id = NA,
            compound.formula = NA,
            compound.adduct = NA,
            compound.rt = NA
          )
        }

        return(res)

      })
      node.anno <- rbindlist(node.annos)
      table(node.anno$type)




      cn.seed.vdata.anned <- cn.seed.vdata%>%
        dplyr::left_join(node.anno,by = "name")

      adducts.score <- adducts %>%
        dplyr::ungroup()%>%
        dplyr::arrange(abs(Mass))%>%
        dplyr::mutate(score = (1:n()),
                      score = 1- score/max(score))%>%
        dplyr::pull(score,name = Adduct)

      cn.seed.vdata.known <- cn.seed.vdata.anned %>%
        dplyr::filter(type %in% "metabolite")%>%
        dplyr::mutate(score.ad = adducts.score[compound.adduct],
                      score.rt =  1- abs(rt - compound.rt)/1e5 ,
                      score.rt = ifelse(score.rt <0 , 0,score.rt),
                      score = score.rt + score.ad) %>%
        dplyr::group_by(compound.id)%>%
        dplyr::arrange( desc(score) )%>%
        dplyr::mutate(temp = seq_len(n()),
                      type = case_when(
                        temp == 1~ type,
                        T~ "adduct"),
                      seed = head(name,n = 1)
        )

      cn.seed.vdata.fi <- cn.seed.vdata.anned%>%
        dplyr::filter(type %in%  c("fragment","isotope"))


      cn.seed.vdata.unknown <- cn.seed.vdata.anned%>%
        dplyr::filter(type %in% "unknown")%>%
        dplyr::mutate(
          type = ifelse(name %in% cn.seed.vdata.fi$seed, "metabolite",type),
          seed = name
        )

      uk.m <- igraph_filter_vertex(cn.seed.ig,
                                   dplyr::pull(cn.seed.vdata.unknown,name))%>%
        get_igraph_membership()
      cn.seed.vdata.unknown <- cn.seed.vdata.unknown%>%
        dplyr::mutate(membership = uk.m[name])%>%
        dplyr::group_by(membership)%>%
        dplyr::arrange(type)%>%
        dplyr::mutate(temp = seq_len(n()),
                      type = case_when(
                        temp == 1 ~ "metabolite",
                        T~ "adduct"),
                      seed = head(name,n = 1)
        )


      cn.seed.vdata3 <- rbind(cn.seed.vdata.known,
                              cn.seed.vdata.unknown,
                              cn.seed.vdata.fi)%>%
        dplyr::select(feature_id = name,
                      TRACE_formula,
                      mz,rt,type,seed,
                      compound_id = compound.id,compound.formula,
                      compound.adduct,compound.rt)
      trace.res <- left_join(cn.seed.vdata3,
                             cpdb[,c("compound_id","name","kegg_id")],
                             by = "compound_id")

    }




    object@advancedAna$TRACE[[pol]] <- trace.res


  }


  return(object)







}




extract_formula_CN <- function(x){
  x <- stringr::str_extract_all(x, "[CN]\\d+") |>
    vapply(paste0, collapse = "", FUN.VALUE = character(1))
  paste0(x, ifelse(grepl("N", x), "", "N0"))
}



