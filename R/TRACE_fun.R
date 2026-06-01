PAVE2 <- function(object){


  if(F){
    object <- readRDS("temp/20251119.rds")
    i.pol = 0
  }


  for (i.pol in 0:1) {

    rt.tol = 10
    ppm = 10

    ### data
    {

      time.start <- Sys.time()
      pol <- ifelse(i.pol==0,"Negative","Positive")
      message_with_time(pol)

      xcms.xcms <- object@xcmsData[[paste0(pol,"MS1")]]
      xcms.net <- get_xcms_feature_connect(xcms.xcms,rt.tol = rt.tol)
      xcms.val <- featureValues(xcms.xcms,missing = 0,value = "maxo")
      xcms.pave.sample <- pData(xcms.xcms)%>%
        dplyr::filter(sample.type %in% c("S12C14N","S12C15N","S13C14N","S13C15N","Blank"))

    }

    ### theoritical mass diff match
    {

      cn.mass.diff <- get_CN_mass_diff_table(N_max = 10)[,type :="CN_label"]
      ad.mass.diff <- get_adduct_mass_diff(unique(polarity(xcms.xcms)))[,type := "adduct"]
      is.mass.diff <- get_iso_mass_diff()[,type := "isotope"]
      fg.mass.diff <- get_fragment_mass_diff()[,type := "fragment"]

      mass.diff.range <- range(cn.mass.diff$mass_diff,
                               ad.mass.diff$mass_diff,
                               is.mass.diff$mass_diff)
      xcms.net <- xcms.net[between.range(mz.diff,  mass.diff.range)]

      cn.match <- match_mz_foverlaps( xcms.net$mz.diff,cn.mass.diff$mass_diff,
                                      ppm.base = xcms.net$mz.mean,ppm = ppm)

      ad.match <- match_mz_foverlaps( xcms.net$mz.diff,ad.mass.diff$mass_diff,
                                      ppm.base = xcms.net$mz.mean,ppm = ppm)

      is.match <- match_mz_foverlaps( xcms.net$mz.diff,is.mass.diff$mass_diff,
                                      ppm.base = xcms.net$mz.mean,ppm = ppm)

      fg.match <- match_mz_foverlaps( xcms.net$mz.diff,fg.mass.diff$mass_diff,
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
        dplyr::mutate(pave_pattern = paste0("C",C_count ,"N",N_count  ))
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
        possible.c.count <- possible.c.count[possible.c.count < c.max&possible.c.count > 0]
        possible.n.count <- possible.n.count[possible.n.count < c.max]
        cn.comb <- expand.grid(C = possible.c.count,
                               N = possible.n.count,
                               p.cor = NA)

        if (!nrow(cn.comb)) return(NULL)
        for (i.cn in 1:nrow(cn.comb)) {
          this.c <- cn.comb$C[i.cn]
          this.n <- cn.comb$N[i.cn]
          all.form <- c(paste0("C0N",this.n,""),paste0("C",this.c,"N0"),paste0("C",this.c,"N",this.n,""))
          all.form <- setdiff(all.form,"C0N0")
          if (!all(all.form %in% x.cn$pave_pattern ) ) next

          #message(x)
          to.id <- x.cn$to[match(all.form,x.cn$pave_pattern)]
          m.detected <- xcms.val[c(x.cn$from[1],to.id),  xcms.pave.sample$sampleNames]
          colnames(m.detected) <- xcms.pave.sample$sample.type
          rownames(m.detected) <- c("C0N0",all.form)
          mean.c0n0 <- mean(m.detected[rownames(m.detected) == "C0N0",
                                       colnames(m.detected) == "S12C14N"])
          m.detected <- m.detected/mean.c0n0
          m.ideal <- get_ideal_CN_ratio(this.c,this.n)%>%t
          m.ideal <- m.ideal[rownames(m.detected),colnames(m.detected)]
          p.cor <- cor(as.vector(m.detected),as.vector(m.ideal))
          cn.comb$p.cor[i.cn] <- p.cor
        }

        p.cor.max <- max(cn.comb$p.cor,na.rm = T)
        #message(p.cor.max)
        if(p.cor.max< 0) return(NULL)
        cn.comb <- cn.comb%>%dplyr::slice_max(p.cor,with_ties = F)
        all.form <- c(paste0("C0N",cn.comb$N,""),paste0("C",cn.comb$C,"N0"),
                      paste0("C",cn.comb$C,"N",cn.comb$N,""))
        all.form <- setdiff(all.form,"C0N0")
        x.cn <- x.cn[match(all.form,x.cn$pave_pattern),]
        x.cn$pave_cor <- p.cor.max
        x.cn$pave_formula <-  paste0("C",cn.comb$C,"N",cn.comb$N,"")
        return(x.cn)

      },BPPARAM = SerialParam(progressbar = T))
      names(cn.net.list.hit) <- names(cn.net.list)
      cn.net.list.hit <- cn.net.list.hit[!sapply(cn.net.list.hit,is.null)]
      cn.net.hit <- data.table::rbindlist(cn.net.list.hit)%>%
        dplyr::filter(pave_cor > 0.75)




      ### timer and temp save
      {
        time.end <- Sys.time()
        time.cost <- difftime(time.end,time.start,units = "secs")%>%as.numeric()
        readr::write_lines(
          paste0("pave2",",",
                 nrow(xcms::featureDefinitions(xcms.xcms)),",",
                 time.cost),
          append = T,
          file = get_dir_expand_from_onedrive("Documents/YLF_Lab/Project/2025.10.10.PAVE/data/pave.CN.count.timer.csv")
        )

        if (T) {

          cn.temp <- data.table::rbindlist(cn.net.list.hit)
          object@advancedAna$PAVE2_temp[[pol]]$CNfinder <- cn.net %>%
            dplyr::mutate(pave_cor = cn.temp$pave_cor[match(ion1,cn.temp$ion1)])

        }

        }
    }


    ### RT and mz error evaluation
    if(T){

      ppm.dyn <- mad(cn.net.hit$mz.ppm) * qnorm(0.99)
      rt.tol.dyn <- mad(cn.net.hit$rt.diff) * qnorm(0.99)

      cn.net.eval <- cn.net%>%
        dplyr::mutate(cn.hit = ion1%in% cn.net.hit$ion1)%>%
        dplyr::arrange(cn.hit)
      data.table::setDT(cn.net.eval)
      cn.net.eval[(cn.hit),mz.ppm]
      #cn.net.eval <- cn.net.eval[1:1000000,]
      ppm.fit <- distinct_norm_from_random_backgroud(cn.net.eval[(cn.hit),mz.ppm],
                                                     cn.net.eval[!(cn.hit),mz.ppm])
      ppm.dyn <- ppm.fit$sd * qnorm(0.99)
      rt.fit <- distinct_norm_from_random_backgroud(cn.net.eval[(cn.hit),rt.diff],
                                                    cn.net.eval[!(cn.hit),rt.diff])
      rt.tol.dyn <- rt.fit$sd * qnorm(0.99)

      object@advancedAna$PAVE2_temp[[pol]][["mz.dyn"]] <- ppm.fit
      object@advancedAna$PAVE2_temp[[pol]][["rt.dyn"]] <- rt.fit

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
      fo <- paste0(object@projectInfo$MSdevFile,".pave.error.pdf")
      fo <- paste0( "C:\\Users\\91879\\OneDrive\\Documents\\YLF_Lab\\Project\\2025.10.10.PAVE\\result/dynamic error/",
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



    ### integration
    {


      xcms.net.candidate <- bind_rows(ad.net,is.net,cn.net.hit,fg.net) %>% data.table::as.data.table()
      xcms.net.candidate <- xcms.net.candidate[
        ,chemform_diff := chemform_simplify(chemform_diff)][,is_CN := (type == "CN_label")
        ][,temp := data.table::fcase(
          is_CN , pave_pattern,
          type == "adduct",paste0(adduct.from,">>",adduct.to),
          default = chemform_diff)
        ][,label :=  paste0(type,": ",temp)
        ]

      ### select unique chemform_diff
      xcms.net.candidate <- xcms.net.candidate[
        ,temp := factor(type, levels = c("CN_label","fragment","adduct","isotope"))
      ][order(temp), .SD[1], by = .(ion1,chemform_diff)]
      data.table::setcolorder(xcms.net.candidate,c("from","to"))
      xcms.net.candidate <- xcms.net.candidate[,eid :=1:.N]

    }


    ### graph of cn seed
    {

      cn.seed.formula <- cn.net.hit %>%
        dplyr::mutate(seed = as.character(from))%>%
        dplyr::distinct(seed, pave_formula)%>%
        dplyr::pull(pave_formula,name = "seed")

      xcms.ig <- igraph::graph_from_data_frame(xcms.net.candidate)
      vda <- vdata(xcms.ig)%>%
        dplyr::mutate(
          color = case_when(
            name %in% cn.net.hit$from ~ "#E64B35",
            T~"#97C2FC"  ),
          pave_formula = cn.seed.formula[name],
          label = paste0(pave_formula,"\n",name))
      vda -> vdata(xcms.ig)


    }


    ### annotate CN seed network edge
    {

      cn.seed <- as.character( cn.net.hit$from)
      cn.seed.ig <- igraph_filter_vertex(xcms.ig,cn.seed)
      .equal.cn.chemform <- function(cn.diff,chemfrom.diff){
        m <- MSCC::chemform_parse(c(cn.diff,chemfrom.diff))
        m <- get_matrix_value_fill_with_NA(m,colnames_vec = c("C","N"))
        m[is.na(m)] <- 0
        m1 <- m[seq_along(cn.diff),]
        m2 <- m[length(cn.diff)+seq_along(cn.diff),]
        eq <- m1[,"C"]==m2[,"C"]&m1[,"N"]==m2[,"N"]
        unname(eq)
      }
      cn.seed.net <- edata(cn.seed.ig)%>%
        dplyr::mutate(
          from.cn = cn.seed.formula[as.character(from)] ,
          to.cn = cn.seed.formula[as.character(to)] ,
          cn.diff = MSCC::chemform_calc(to.cn,from.cn ,"-",return = "chemform"),
          cn.equal = cn.diff == "",
          chemform.equal = .equal.cn.chemform(cn.diff,chemform_diff)
        )
      cn.seed.net <- cn.seed.net%>%
        dplyr::mutate(
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

          )
        )%>%data.table()



      if(F){

        for (i.csi in 1:nrow(cn.seed.net)) {

          i.csi.from <- cn.seed.net$from[i.csi]
          i.csi.to <- cn.seed.net$to[i.csi]
          from.peaks <- cn.net.list.hit[[as.character(i.csi.from)]]
          to.peaks <- cn.net.list.hit[[as.character(i.csi.to)]]

          from.matrix <- xcms.val[c(from.peaks$from[1],from.peaks$to),xcms.pave.sample$sampleNames]
          dimnames(from.matrix) <- list(c("C0N0",from.peaks$pave_pattern),xcms.pave.sample$sample.type)
          to.matrix <- xcms.val[c(to.peaks$from[1],to.peaks$to),xcms.pave.sample$sampleNames]
          dimnames(to.matrix) <- list(c("C0N0",to.peaks$pave_pattern),xcms.pave.sample$sample.type)
          C0N0.mean <- mean(from.matrix["C0N0",colnames(from.matrix) == "S12C14N"])
          from.matrix <- from.matrix/C0N0.mean
          to.matrix <- to.matrix/C0N0.mean
          #cor(as.vector(from.matrix),as.vector(to.matrix))
          if (F) {

            #print(rownames(to.matrix))
            #print(rownames(from.matrix))
            #rowMeans(to.matrix/from.matrix,na.rm = T)%>%print()

            m <- rbind(from.matrix,to.matrix)
            if (nrow(m)!=8) {
              next
            }
            rownames(m) <- c("C0N0","C0Ny","CxN0","CxNy",
                             "C1N0","C1Ny","Cx-1N0","Cx-1Ny")
            hm <- Heatmap(
              m,name = "Relative\nAbundance",
              column_title = paste0("x = ",get_formula_ele_count(rownames(from.matrix)[4],"C")," , ",
                                    "y = ",get_formula_ele_count(rownames(from.matrix)[4],"N")  ),
              cluster_rows = F,cluster_columns = F,
              col = colramp(),
              column_names_rot = -30,
              row_names_side = "left")
            open_plot_win(hm,8,5)
          }

        }
      }


    }



    ### graph annotate CN seed, using Loop to determine
    if(F){

      message_with_time("Graph annotate CN seed")
      cn.seed.ig <- igraph_filter_vertex(xcms.ig,cn.seed)
      cn.seed.ig <- igraph_remove_edge(cn.seed.ig,
                                       which(E(cn.seed.ig)$eid %in% cn.seed.net[new.type=="false",eid]))

      cn.seed.vdata <- vdata(cn.seed.ig)%>%
        dplyr::mutate(node.group = get_igraph_membership(cn.seed.ig))
      cn.seed.split <- split(cn.seed.vdata$name,cn.seed.vdata$node.group)
      table(lengths(cn.seed.split))
      cn.seed.annotation <- list()
      pb <- get_progress_bar(total_iterations =  length(cn.seed.split))
      for (i.css in seq_along(cn.seed.split) ) {

        pb$tick()
        #i.css <- which(lengths(cn.seed.split) == 6 )[2]
        i.seed <- cn.seed.split[[i.css]]
        #message_with_time(i.css,"; ",length(i.seed))
        i.cn.seed.ig <- igraph_filter_vertex(cn.seed.ig, i.seed)
        #vis_pave_igraph(i.cn.seed.ig)

        if (length(i.seed) == 1) {

          cn.seed.annotation[[i.css]] <- data.frame(
            name = i.seed,
            pave_MS_form = "Undefined",
            pave_formula = cn.seed.formula[i.seed],
            pave_seed = paste0("CN_Seed_",i.seed)

          )
          next
        }

        i.cn.seed.ig.ring <- igraph::simple_cycles(
          i.cn.seed.ig,mode = "all",min = 3,max = 4)

        if (length(i.cn.seed.ig.ring$vertices)>0) {
          ### ring test
          {

            eloop <- list()
            ring.node.forms <- list()
            for(i.ring in seq_len(length(i.cn.seed.ig.ring$vertices)) ){

              i.ep <-  i.cn.seed.ig.ring$edges[[i.ring]]
              i.dir <- get_path_direction(i.cn.seed.ig,
                                          i.cn.seed.ig.ring$vertices[[i.ring]],
                                          i.ep)
              i.cd <- i.cn.seed.ig.ring$edges[[i.ring]]$chemform_diff
              i.cd <- MSCC::chemform_multi(i.cd,i.dir,return = "chemform")
              #i.cd.cum <- MSCC:::chemform_sum(i.cd)
              i.cd.temp <- sapply(seq_along(i.cd),function(x){
                MSCC:::chemform_sum(i.cd[1:x])
              })
              i.cd.temp <- chemform_remove_iso(i.cd.temp)
              i.cd.temp <- chemform_simplify(i.cd.temp)
              i.cd.exist <- (i.cd.temp%in%  c(fg.mass.diff$chemform_diff,
                                              ad.mass.diff$chemform_diff))

              if (all(i.cd.exist)) {

                i.ring.ig  <- igraph_filter_edge(i.cn.seed.ig,
                                                 which(E(i.cn.seed.ig)$eid %in%  i.ep$eid))
                i.ring.vform <- get_pave_ig_vertex_form(i.ring.ig)
                if (any(lengths(i.ring.vform) > 1)) next
                i.ring.vform <- unlist(i.ring.vform)
                #print(i.ring.vform)
                ring.node.forms[[i.ring]] <- i.ring.vform
                eloop[[i.ring]] <- i.ep$eid
              }
            }



          }

          ### determine adduct frag iso form
          if(length(ring.node.forms)>0){


            ring.node.form <- do.call(bind_rows, ring.node.forms)
            rnfg <- ring.node.form.group(ring.node.form)
            rnfc <- ring.node.form%>%
              dplyr::mutate(id= 1:n())%>%
              tidyr::pivot_longer(!id)%>%
              dplyr::mutate(cl = rnfg[id])%>%
              dplyr::filter(!is.na(value))%>%
              dplyr::group_by(name)%>%
              dplyr::mutate(cl.freq = names(which.max(table(cl))))%>%
              dplyr::filter(cl == cl.freq)%>%
              dplyr::ungroup()%>%
              dplyr::distinct(name,value,cl)%>%
              dplyr::group_by(cl)%>%
              dplyr::mutate(cl = paste0("CN_Seed_",min(name)))
            vc <- dplyr::pull(rnfc,cl,name= name)
            vf <- dplyr::pull(rnfc,value,name= name)
            eloop <- unique(unlist(eloop))


            ### annotate
            {

              cn.seed.annotation[[i.css]] <- data.frame(
                name = i.seed
              ) %>%
                dplyr::mutate(
                  pave_formula = cn.seed.formula[i.seed],
                  pave_MS_form = case_when(
                    name %in%names(vf)~ vf[name],
                    T~"Unknow MS form"),
                  pave_seed = case_when(
                    name %in%names(vc)~ vc[name],
                    T~ NA)
                )%>%dplyr::select(name,pave_formula,pave_MS_form,pave_seed)

              next

            }

          }

        }

        if (T){

          if (all(E(i.cn.seed.ig)$type=="fragment")) {
            i.seed.seed  <- edata(i.cn.seed.ig)$to[1]
          }else{
            i.cn.seed.ig.remove.fg <- igraph_filter_edge(
              i.cn.seed.ig, which(E(i.cn.seed.ig)$type!="fragment"))
            i.seed.seed <- as.character(min(as.numeric(names(V(i.cn.seed.ig.remove.fg)))))
          }

          i.seed.form <- get_pave_ig_vertex_form(i.cn.seed.ig)%>%
            sapply(`[`,1)
          #i.seed.form <- unlist(i.seed.form)
          if (length(i.seed.form) > length(i.seed)) {
            break
          }

          cn.seed.annotation[[i.css]] <- data.frame(
            name = i.seed,
            pave_MS_form = i.seed.form[i.seed],
            pave_formula = cn.seed.formula[i.seed],
            pave_seed = paste0("CN_Seed_",i.seed.seed)
          )
          next
        }







      }

      pave.cor <- dplyr::pull(.data = cn.net.hit,pave_cor,name = from)
      cn.seed.annotation.df <- do.call(bind_rows,cn.seed.annotation)%>%
        dplyr::mutate(
          pave_annotation = case_when(
            grepl(pattern = "^\\[",pave_MS_form)~"isotope",
            pave_seed == paste0("CN_Seed_",name)~ "CN_metabolite",
            is.na(pave_seed) ~"Unknow_adduct",
            !grepl(pattern = ";$",pave_MS_form)~"fragment",
            T~"adduct"
          ),
          pave_cor = pave.cor[name],
          pave_pattern = "C0N0",
          pave_cn_seed = name

        )%>%data.table::setDT()

      cn.exp <- cn.net.hit[
        from  %in% cn.seed.annotation.df$name][
          , name := as.character(to) ][
            ,pave_cn_seed := as.character(from)][
              , .(name,pave_formula   ,pave_pattern  ,pave_cor ,pave_cn_seed)
            ][cn.seed.annotation.df[,.(pave_seed,pave_MS_form,pave_annotation,pave_cn_seed)],
              on = "pave_cn_seed"]

      cn.peaks.annotation.df <- bind_rows(cn.seed.annotation.df,cn.exp)
      table(cn.peaks.annotation.df$pave_annotation)

    }


    ### graph annotate CN seed
    if(F){

      message_with_time("Graph annotate CN seed")
      cn.seed.ig <- igraph_filter_vertex(xcms.ig,cn.seed)
      cn.seed.ig <- igraph_remove_edge(cn.seed.ig,
                                       which(E(cn.seed.ig)$eid %in% cn.seed.net[new.type=="false",eid]))

      cn.seed.vdata <- vdata(cn.seed.ig)%>%
        dplyr::mutate(node.group = get_igraph_membership(cn.seed.ig))
      cn.seed.split <- split(cn.seed.vdata$name,cn.seed.vdata$node.group)
      table(lengths(cn.seed.split))
      cn.seed.annotation <- list()
      for (i.css in seq_along(cn.seed.split) ) {

        #i.css <- which(lengths(cn.seed.split) == 12 )[1]
        i.seed <- cn.seed.split[[i.css]]
        #message_with_time(i.css,"; ",length(i.seed))
        i.cn.seed.ig <- igraph_filter_vertex(cn.seed.ig, i.seed)
        vis_pave_igraph(i.cn.seed.ig)
        vis_pave_igraph(pave_igraph_contract(i.cn.seed.ig))

        if (length(i.seed) == 1) {

          cn.seed.annotation[[i.css]] <- data.frame(
            name = i.seed,
            pave_MS_form = "Undefined",
            pave_formula = cn.seed.formula[i.seed],
            pave_seed = paste0("CN_Seed_",i.seed)

          )
          next
        }
        if (length(i.seed) > 1) {

          i.cn.seed.ig.ct <- pave_igraph_contract(i.cn.seed.ig)
          eda <- edata(i.cn.seed.ig.ct) %>%
            dplyr::filter(!from == to)

          if (nrow(eda) == 0) {

            if (all(E(i.cn.seed.ig)$type=="fragment")) {
              i.seed.seed  <- edata(i.cn.seed.ig)$to[1]
            }else{
              i.cn.seed.ig.remove.fg <- igraph_filter_edge(
                i.cn.seed.ig, which(E(i.cn.seed.ig)$type!="fragment"))
              i.seed.seed <- as.character(min(as.numeric(names(V(i.cn.seed.ig.remove.fg)))))
            }

            i.seed.form <- get_pave_ig_vertex_form(i.cn.seed.ig)%>%
              sapply(`[`,1)
            #i.seed.form <- unlist(i.seed.form)
            if (length(i.seed.form) > length(i.seed)) {
              break
            }

            cn.seed.annotation[[i.css]] <- data.frame(
              name = i.seed,
              pave_MS_form = i.seed.form[i.seed],
              pave_formula = cn.seed.formula[i.seed],
              pave_seed = paste0("CN_Seed_",i.seed.seed)
            )
            next
          }

          v.anno <- data.table(
            name = c(eda$from,eda$to),
            adduct = c(eda$adduct.from,eda$adduct.to),
            eid = c(eda$eid,eda$eid),
            error = c(abs(eda$rt.diff),abs(eda$rt.diff))
          )%>%
            dplyr::group_by(name)%>%
            dplyr::mutate(
              degree =n() )%>%
            dplyr::group_by(name,adduct)%>%
            dplyr::mutate(freq = n(),
                          ratio = freq/ degree)%>%
            dplyr::ungroup()%>%
            dplyr::mutate(
              ratio = case_when(degree == 1~ 0.5,
                                T~ratio)
            )

          eid.filt <- names(which(sort(mean_f(v.anno$ratio,v.anno$eid)) > 0.5))
          igraph_filter_edge(i.cn.seed.ig.ct, which(E(i.cn.seed.ig.ct)$eid %in%eid.filt))%>%
            vis_pave_igraph()


          eid.reomve <- names(which(sort(mean_f(v.anno$ratio,v.anno$eid)) <= 0.5))
          igraph_remove_edge(i.cn.seed.ig, which(E(i.cn.seed.ig.ct)$eid %in%eid.reomve))%>%
            vis_pave_igraph()

        }






      }

      cn.seed.annotation.df <- do.call(bind_rows,cn.seed.annotation)%>%
        dplyr::mutate(
          pave_annotation = case_when(
            grepl(pattern = "^\\[",pave_MS_form)~"isotope",
            pave_seed == paste0("CN_Seed_",name)~ "CN_metabolite",
            is.na(pave_seed) ~"Unknow_adduct",
            !grepl(pattern = ";$",pave_MS_form)~"fragment",
            T~"adduct"
          )
        )

      rownames(cn.seed.annotation.df) <- cn.seed.annotation.df$name

      cn.exp <- cn.seed.annotation.df[as.character(cn.net.hit$from),]%>%
        dplyr::mutate(name = as.character(cn.net.hit$to))
      cn.seed.annotation.df <- rbind(cn.seed.annotation.df,cn.exp)
      table(cn.seed.annotation.df$pave_annotation)
    }


    ### formula assign
    {

      cpdb_path <- object@projectInfo$CompoundDB_path
      cpdb <- CompoundDb::CompDb(cpdb_path)
      xcms.fdf <- xcms::featureDefinitions(xcms.xcms)


    }


    object@advancedAna$PAVE2[[pol]] <- cn.peaks.annotation.df


  }


  return(object)







}


distinct_norm_from_random_backgroud <- function(
    x_norm,x_random){


  bg_kde <- density(x_random, n = 1e5)

  f_bg <- function(x) {
    approx(bg_kde$x, bg_kde$y, xout = x,
           rule = 2, ties = mean)$y
  }

  fit_bg_norm_mixture <- function(
    x,
    f_bg,
    max_iter = 1e3,
    tol = 1e-6
  ) {
    n <- length(x)

    pi  <- 0.2
    mu  <- mean(x)
    sd  <- sd(x)

    loglik_old <- -Inf

    for (iter in seq_len(max_iter)) {

      bg_d <- f_bg(x)
      bg_d[bg_d <= 0] <- min(bg_d[bg_d > 0]) * 1e-3

      norm_d <- dnorm(x, mu, sd)

      w <- pi * norm_d / ((1 - pi) * bg_d + pi * norm_d)

      pi <- mean(w)
      mu <- sum(w * x) / sum(w)
      sd <- sqrt(sum(w * (x - mu)^2) / sum(w))

      loglik <- sum(log((1 - pi) * bg_d + pi * norm_d))
      if (abs(loglik - loglik_old) < tol) break
      loglik_old <- loglik
    }

    list(
      pi = pi,
      mu = mu,
      sd = sd,
      posterior = w,
      iter = iter
    )
  }

  fit <- fit_bg_norm_mixture(x_norm, f_bg)
  return(fit)

}
