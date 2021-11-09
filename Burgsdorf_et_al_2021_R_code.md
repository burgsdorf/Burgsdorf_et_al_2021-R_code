R code used to build figures for Burgsdorf et al., 2021
================

## Lineage-specific energy and carbon metabolism of sponge symbionts and contributions to the host carbon pool

### Final figures were graphically edited using Inkscape (<a href="http://inkscape.org" class="uri">http://inkscape.org</a>)

## Figure 2

Data used to build the final figure construction.

    load("D:/Dropbox/(Bio)Informatic_protocols/R/!!!_Markdown_pipelines/2020_50_genomes_graphics_UTF-8.RData")
    library(ggplot2)

    ## Warning: package 'ggplot2' was built under R version 3.6.3

    library(superheat)

    ## Warning: package 'superheat' was built under R version 3.6.3

    library(plyr)

    ## Warning: package 'plyr' was built under R version 3.6.3

    library(ggpubr)

    ## Warning: package 'ggpubr' was built under R version 3.6.3

    ## 
    ## Attaching package: 'ggpubr'

    ## The following object is masked from 'package:plyr':
    ## 
    ##     mutate

    library(reshape2)

    ## Warning: package 'reshape2' was built under R version 3.6.3

    # set the text colors 
    # identify all scaled values that fall below (white)
    LifeStyle_RA2.col2 <- LifeStyle_RA2 > 35
    # set all values that satisfy the condition to "white"
    LifeStyle_RA2.col2 <- gsub("TRUE", "black", LifeStyle_RA2.col2)
    # set all values that do not satisfy the condition to "black"
    LifeStyle_RA2.col2 <- gsub("FALSE", "white", LifeStyle_RA2.col2)
    # convert to matrix
    LifeStyle_RA2.col2 <- matrix(LifeStyle_RA2.col2, ncol = ncol(LifeStyle_RA2))

    superheat(LifeStyle_RA2,
              X.text = round(as.matrix(LifeStyle_No_genes2), 2),
              X.text.size = 4.5,
              X.text.col = LifeStyle_RA2.col2,
              left.label.size = 0.4,
              smooth.heat = TRUE,
              left.label.text.size = 5,
              left.label.text.alignment = "right",
              left.label.col = "white",
              bottom.label.col = 'white',
              bottom.label.size = 0.5,
              bottom.label.text.size = 5,
              bottom.label.text.angle = 45,
              bottom.label.text.alignment = "right",
              pretty.order.cols = FALSE,
              pretty.order.rows = TRUE,
              #heat.pal = viridis::inferno(100),
              clustering.method = "hierarchical",
              dist.method = "binary",
              bottom.label.text.col =c("darkmagenta","darkmagenta","chartreuse4","chartreuse4","black","black","black","firebrick2","firebrick2"))

    ## [1]   0  20  50  80 100

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](Burgsdorf_et_al_2021_R_code_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Figure 3

### Figure 3A

    COG1529_taxa_bbplot <- ggplot(COG1529_taxa, aes(x=Taxonomy, y=COG1529, fill=Taxonomy)) +
      ylab("No of COG1529 per genome") +
      xlab("Taxonomy (Phylum/Class)") + 
      geom_boxplot() + 
      coord_flip() + 
      theme_classic()  +
      scale_fill_manual(values=group_pasteur_21.colors) +
      theme(legend.position = "none", text = element_text(size=14), axis.text=element_text(size=14)) +
      scale_y_continuous(breaks = seq(0, 60, 10))

    COG1529_taxa_bbplot

![](Burgsdorf_et_al_2021_R_code_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Figure 3B

    # set the text colors 
    # identify all scaled values that fall below (white)
    COG1529_taxa_RA.col2 <- COG1529_taxa_RA > 1
    # set all values that satisfy the condition to "white"
    COG1529_taxa_RA.col2 <- gsub("TRUE", "black", COG1529_taxa_RA.col2)
    # set all values that do not satisfy the condition to "black"
    COG1529_taxa_RA.col2 <- gsub("FALSE", "white", COG1529_taxa_RA.col2)
    # convert to matrix
    COG1529_taxa_RA.col2 <- matrix(COG1529_taxa_RA.col2, ncol = ncol(COG1529_taxa_RA))

    superheat(COG1529_taxa_RA,
              X.text = round(as.matrix(COG1529_taxa_No_genes), 2),
              X.text.size = 5,
              X.text.col = COG1529_taxa_RA.col2,
              left.label.size = 0.8,
              left.label.text.size = 5,
              left.label.text.alignment = "right",
              left.label.col = "white",
              bottom.label.col = 'white',
              bottom.label.size = 0.3,
              bottom.label.text.size = 5,
              bottom.label.text.angle = 90,
              bottom.label.text.alignment = "right",
              pretty.order.cols = TRUE,
              pretty.order.rows = TRUE,
              yr = COG1529_taxa_RA$K03520,
              yr.axis.name = "K03520 (%)",
              yr.plot.type = "scatterline",
              yr.line.col = "tomato3",
              yr.obs.col = rep("orange", nrow(COG1529_taxa_RA)),
              yr.point.size = 4,
              # order the rows by mpg
              order.rows = order(COG1529_taxa_RA$K03520)
              )

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

    ## Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique
    ## 'x' values

![](Burgsdorf_et_al_2021_R_code_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Figure 5

### Preliminary data before binary transformation. An example of *cox* genes.

    cox_ID_THRH_melted_g1 <- ggplot(cox_ID_THRH_melted, aes(x = `Final annotation`, y = value_sqrt_sqrt, colour = `Phylum / Class`)) +
      geom_jitter(alpha=0.25, show.legend = FALSE) +
      scale_color_manual(values=group.colors) +
      xlab("") +
      ylab("Expression (root 4th transformed)") +
      stat_summary(aes(colour = `Phylum / Class`), fun = mean, geom = "point", size = 5, alpha = 1, position = pd2) +
      #scale_y_continuous(breaks = seq(0, 6, 1)) +
      theme_bw() +
      theme(legend.position = "bottom", legend.direction = "horizontal", text = element_text(size=15))
      
    cox_ID_THRH_melted_g1

![](Burgsdorf_et_al_2021_R_code_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

    cox_ID_THRH_short2_agg_sep_melted <- melt(cox_ID_THRH_short2_agg_sep[,-2],id=c("Final_annotation","Phylum/Class"))
    cox_ID_THRH_short2_agg_sep_melted$binary <- ifelse(cox_ID_THRH_short2_agg_sep_melted$value>0, 1, 0)

    cox_ID_THRH_melted_gF2 <- ggplot(cox_ID_THRH_short2_agg_sep_melted, aes(x = `Final_annotation`, y = binary, colour = `Phylum/Class`)) +
      geom_point(position=position_jitter(h=0.15, w=0.45), shape = 20, size = 3, alpha=0.25) +
      scale_color_manual(values=group.colors) +
      stat_summary(data = cox_ID_THRH_short2_NA_sep_melted, aes(x = `Final_annotation`, y = binary), fun = mean, geom = "point", size = 5, alpha = 1) +
      xlab("") +
      ylab("") +
      stat_summary(aes(colour = `Phylum/Class`), fun = mean, geom = "point", size = 5, alpha = 1, position = pd2) +
      stat_summary(aes(group = `Phylum/Class`), fun = mean, geom = "line", position = pd2) + 
      scale_y_continuous(breaks = seq(0, 1, 0.5)) +
      theme_bw() +
      theme(legend.position = "none", legend.direction = "horizontal", text = element_text(size=15))

    anaplerotic_4_genes_ID_THRH_melted_gF3 <- ggplot(anaplerotic_4_genes_ID_THRH_short2_agg_sep_melted, aes(x = `Final_annotation`, y = binary, colour = `Phylum/Class`)) +
      geom_point(position=position_jitter(h=0.15, w=0.45), shape = 20, size = 3, alpha=0.25) +
      scale_color_manual(values=group.colors) +
      xlab("") +
      ylab("Expression (binary)") +
      stat_summary(aes(colour = `Phylum/Class`), fun = mean, geom = "point", size = 5, alpha = 1, position = pd2) +
      stat_summary(data = anaplerotic_4_genes_ID_THRH_short2_NA_sep_melted, aes(x = `Final_annotation`, y = binary), fun = mean, geom = "point", size = 5, alpha = 1) +
      scale_y_continuous(breaks = seq(0, 1, 0.5)) +
      theme_bw() +
      theme(legend.position = "none", legend.direction = "horizontal", text = element_text(size=15))

    fixation_3_ID_THRH_melted_gF3 <- ggplot(fixation_3_ID_THRH_short2_agg_sep_melted, aes(x = `Final_annotation`, y = binary, colour = `Phylum/Class`)) +
      geom_point(position=position_jitter(h=0.15, w=0.45), shape = 20, size = 3, alpha=0.25) +
      scale_color_manual(values=group.colors) +
      xlab("") +
      ylab("") +
      stat_summary(aes(colour = `Phylum/Class`), fun = mean, geom = "point", size = 5, alpha = 1, position = pd2) +
      stat_summary(data = fixation_3_ID_THRH_short2_NA_sep_melted, aes(x = `Final_annotation`, y = binary), fun = mean, geom = "point", size = 5, alpha = 1) +
      scale_y_continuous(breaks = seq(0, 1, 0.5)) +
      theme_bw() +
      theme(legend.position = "none", legend.direction = "horizontal", text = element_text(size=15))

    ggarrange(cox_ID_THRH_melted_gF2, anaplerotic_4_genes_ID_THRH_melted_gF3, fixation_3_ID_THRH_melted_gF3,
              labels = c("A", "B", "C"),
              ncol = 1, nrow = 3)

![](Burgsdorf_et_al_2021_R_code_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Figure 6

    load("D:/Dropbox/(Bio)Informatic_protocols/R/!!!_Markdown_pipelines/Radioactive_sponge_2020.RData")

### An example of mean and sd calculation based on “status” (dark or light) and core depths (mm) of Petrosia ficiformis.

    Pf_Exp2_ddply <- ddply(Pf_Exp2, c("status", "mm"), summarise,
                   N    = length(fixed),
                   mean = mean(fixed),
                   sd   = sd(fixed),
                   se   = sd / sqrt(N)
    )

### Creating plot with sd as error bars.

    Ts_Fig12_graph <- ggplot()+
      geom_line(data=Ts_Fig12, aes(x=mm, y=fixed, group=status, color=status), linetype = "dashed", size=1) +
      geom_point(data=Ts_Fig12, aes(x=mm, y=fixed, group=status, color=status), size = 4.5) +
      geom_errorbar(data=Ts_Fig12, aes(x=mm, group=status, color=status, ymin=fixed-SD, ymax=fixed+SD), width=0.15) +
      theme(#legend.position="bottom",
            legend.position = c(0.95, 0.75),
            legend.justification = c(1, 0.1),
            legend.title = element_blank(),
            legend.text=element_text(size=16),
            legend.key=element_blank(),
            plot.margin = unit(c(1,1,1,1), "cm"),
            axis.title.x=element_text(colour="black", size = 17),
            axis.title.y=element_text(colour="black", size = 17),
            panel.background = element_blank(),
            axis.line.x = element_line(size = 0.5, colour = "gray38"),
            axis.line.y = element_line(size = 0.5, colour = "gray38"),
            axis.text.x=element_text(colour="black", size = 17),
            axis.text.y=element_text(colour="black", size = 17)) +
      scale_y_continuous(expand = c(0,0), breaks = seq(0, 50, 5), limits = c(-0.4, 52)) +
      scale_x_continuous(expand = c(0,0), breaks = seq(0, 14, 2), limits = c(0, 15)) +
      scale_color_manual(labels = c("+CaSs Dark", "+CaSs Light"), values = c("royalblue3", "orange2")) +
      xlab("") + ylab("Ci fixed (µg/g)")


    Pf_Exp2_ddply_graph <- ggplot()+
      geom_line(data=Pf_Exp2_ddply, aes(x=mm, y=mean, group=status, color=status), linetype = "dashed", size=1) +
      geom_point(data=Pf_Exp2_ddply, aes(x=mm, y=mean, group=status, color=status), size = 4.5) +
      geom_errorbar(data=Pf_Exp2_ddply, aes(x=mm, group=status, color=status, ymin=mean-sd, ymax=mean+sd), width=0.15) +
      theme(#legend.position="bottom",
            legend.position = c(0.95, 0.75),
            legend.justification = c(1, 0.1),
            legend.title = element_blank(),
            legend.text=element_text(size=16),
            legend.key=element_blank(),
            plot.margin = unit(c(1,1,1,1), "cm"),
            axis.title.x=element_text(colour="black", size = 17),
            axis.title.y=element_text(colour="black", size = 17),
            panel.background = element_blank(),
            axis.line.x = element_line(size = 0.5, colour = "gray38"),
            axis.line.y = element_line(size = 0.5, colour = "gray38"),
            axis.text.x=element_text(colour="black", size = 17),
            axis.text.y=element_text(colour="black", size = 17)) +
      scale_y_continuous(expand = c(0,0), breaks = seq(0, 50, 5), limits = c(-0.4, 52)) +
      scale_x_continuous(expand = c(0,0), breaks = seq(0, 6, 2), limits = c(0, 7)) +
      scale_color_manual(labels = c("+CaSf Dark","-CaSf Dark","+CaSf Light","-CaSf Light"), values = c("royalblue3","royalblue4","orange2","orange3")) +
      xlab("") + ylab("")

    Ts_Fig13_graph <- ggplot()+
      geom_line(data=Ts_Fig13, aes(x=mm, y=fixed, group=status, color=status), linetype = "dashed", size=1) +
      geom_point(data=Ts_Fig13, aes(x=mm, y=fixed, group=status, color=status), size = 4.5) +
      geom_errorbar(data=Ts_Fig13, aes(x=mm, group=status, color=status, ymin=fixed-SD, ymax=fixed+SD), width=0.15) +
      theme(#legend.position="bottom",
            legend.position = c(0.95, 0.75),
            legend.justification = c(1, 0.1),
            legend.title = element_blank(),
            legend.text=element_text(size=16),
            legend.key=element_blank(),
            plot.margin = unit(c(1,1,1,1), "cm"),
            axis.title.x=element_text(colour="black", size = 17),
            axis.title.y=element_text(colour="black", size = 17),
            panel.background = element_blank(),
            axis.line.x = element_line(size = 0.5, colour = "gray38"),
            axis.line.y = element_line(size = 0.5, colour = "gray38"),
            axis.text.x=element_text(colour="black", size = 17),
            axis.text.y=element_text(colour="black", size = 17)) +
      scale_y_continuous(expand = c(0,0), breaks = seq(0, 25, 5), limits = c(-0.4, 27.5)) +
      scale_x_continuous(expand = c(0,0), breaks = seq(0, 14, 2), limits = c(0, 15)) +
      scale_color_manual(labels = c("+CaSs Dark", "+CaSs Light"), values = c("royalblue3", "orange2")) +
      xlab("Sponge section (mm)") + ylab("Ci fixed (µg/g)")

    Pf_Exp3_ddply_graph <- ggplot()+
      geom_line(data=Pf_Exp3_ddply, aes(x=mm, y=mean, group=status, color=status), linetype = "dashed", size=1) +
      geom_point(data=Pf_Exp3_ddply, aes(x=mm, y=mean, group=status, color=status), size = 4.5) +
      geom_errorbar(data=Pf_Exp3_ddply, aes(x=mm, group=status, color=status, ymin=mean-sd, ymax=mean+sd), width=0.15) +
      theme(#legend.position="bottom",
            legend.position = c(0.95, 0.75),
            legend.justification = c(1, 0.1),
            legend.title = element_blank(),
            legend.text=element_text(size=16),
            legend.key=element_blank(),
            plot.margin = unit(c(1,1,1,1), "cm"),
            axis.title.x=element_text(colour="black", size = 17),
            axis.title.y=element_text(colour="black", size = 17),
            panel.background = element_blank(),
            axis.line.x = element_line(size = 0.5, colour = "gray38"),
            axis.line.y = element_line(size = 0.5, colour = "gray38"),
            axis.text.x=element_text(colour="black", size = 17),
            axis.text.y=element_text(colour="black", size = 17)) +
      scale_y_continuous(expand = c(0,0), breaks = seq(0, 25, 5), limits = c(-0.4, 27.5)) +
      scale_x_continuous(expand = c(0,0), breaks = seq(0, 6, 2), limits = c(0, 7)) +
      scale_color_manual(labels = c("+CaSf Dark","-CaSf Dark","+CaSf Light","-CaSf Light"), values = c("royalblue3","royalblue4","orange2","orange3")) +
      xlab("Sponge section (mm)") + ylab("")


    # Combining plots
    ggarrange(Ts_Fig12_graph, Pf_Exp2_ddply_graph, Ts_Fig13_graph, Pf_Exp3_ddply_graph, ncol = 2, nrow = 2, labels = c("A","B","C","D"))

![](Burgsdorf_et_al_2021_R_code_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->
