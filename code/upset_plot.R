upset_R <- function(filepath2) {
    filepath <- list.files(filepath2, pattern = ".txt", full.names = TRUE)
    filename <- gsub(".txt", "", basename(filepath))
    cancername <- unique(gsub("_vs.*.", "", filename))

    # filename2 <- gsub(".*._","",filename)

    for (i in 1:length(filepath)) {
        PTmarkers <- read.table(filepath[i], header = T, sep = "\t")
        PTmarkers <- PTmarkers[PTmarkers$pct.1 > 0.05 & PTmarkers$pct.2 > 0.05, ]
        highmarkers <- rownames(PTmarkers[PTmarkers$p_val_adj < 0.05 & PTmarkers$avg_log2FC > 0.5, ]) ## Putting the cutoff < 0.05 and avg_log2FC > 1
        assign(paste(filename[i], "_sig_high", sep = ""), highmarkers)
        lowmarkers <- rownames(PTmarkers[PTmarkers$p_val_adj < 0.05 & PTmarkers$avg_log2FC < -0.5, ])
        assign(paste(filename[i], "_sig_low", sep = ""), lowmarkers)
    }

    filename2 <- grep(cancername, ls(pattern = "_sig_high"), value = TRUE)

    H <- list(
        "BC_BCC_high" = get(filename2[1]),
        "BC_CM_high" = get(filename2[2]),
        "BC_CRC_high" = get(filename2[3]),
        "BC_EA_high" = get(filename2[4]),
        "BC_HCC_high" = get(filename2[5]),
        "BC_ICC_high" = get(filename2[6]),
        "BC_NSCLC_high" = get(filename2[7]),
        "BC_PDAC_high" = get(filename2[8]),
        "BC_RCC_high" = get(filename2[9]),
        "BC_SCC_high" = get(filename2[10]),
        "BC_UM_high" = get(filename2[11])
    )

    # create customised venn diagram
    # library(ggvenn)
    # a <- ggvenn(H, show_elements = FALSE, stroke_color = "Black", stroke_linetype = "solid")
    # savedir <- "/diazlab/data3/.abhinav/projects/SHH/snRNA/removed_samples_BO/scdifferential/"
    # dir.create(paste(savedir, "Table/", sep = ""), showWarnings = FALSE)
    # pdf(paste(savedir, "Table/MAST_recurrent_high_venny.pdf", sep = ""))
    # a
    # dev.off()

    ### Making an upset plot
    n <- max(
        length(get(filename2[1])), length(get(filename2[2])),
        length(get(filename2[3])), length(get(filename2[4])),
        length(get(filename2[5])), length(get(filename2[6])),
        length(get(filename2[7])), length(get(filename2[8])),
        length(get(filename2[9])), length(get(filename2[10])),
        length(get(filename2[11]))
    )

    length(all_BC_vs_BCC_sig_high) <- n
    length(all_BC_vs_CM_sig_high) <- n
    length(all_BC_vs_CRC_sig_high) <- n
    length(all_BC_vs_EA_sig_high) <- n
    length(all_BC_vs_HCC_sig_high) <- n
    length(all_BC_vs_ICC_sig_high) <- n
    length(all_BC_vs_NSCLC_sig_high) <- n
    length(all_BC_vs_PDAC_sig_high) <- n
    length(all_BC_vs_RCC_sig_high) <- n
    length(all_BC_vs_SCC_sig_high) <- n
    length(all_BC_vs_UM_sig_high) <- n

    df <- as.data.frame(cbind(
        (all_BC_vs_BCC_sig_high), (all_BC_vs_CM_sig_high),
        (all_BC_vs_CRC_sig_high), (all_BC_vs_EA_sig_high),
        (all_BC_vs_HCC_sig_high), (all_BC_vs_ICC_sig_high),
        (all_BC_vs_NSCLC_sig_high), (all_BC_vs_PDAC_sig_high),
        (all_BC_vs_RCC_sig_high), (all_BC_vs_SCC_sig_high),
        (all_BC_vs_UM_sig_high)
    ))

    library(UpSetR)
    fld <- fromList(as.list(df)) ## this is an intersect function

    colnames(fld) <- c(
        ("all_BC_vs_BCC_sig_high"), ("all_BC_vs_CM_sig_high"),
        ("all_BC_vs_CRC_sig_high"), ("all_BC_vs_EA_sig_high"),
        ("all_BC_vs_HCC_sig_high"), ("all_BC_vs_ICC_sig_high"),
        ("all_BC_vs_NSCLC_sig_high"), ("all_BC_vs_PDAC_sig_high"),
        ("all_BC_vs_RCC_sig_high"), ("all_BC_vs_SCC_sig_high"),
        ("all_BC_vs_UM_sig_high")
    )

    dir.create(paste0(filepath2, "/Table"), showWarnings = FALSE, recursive = TRUE)
    pdf(paste(filepath2, "/Table/", cancername, "_all_upset_plot.pdf", sep = ""), width = 20, height = 14)
    upset(fld, nsets = 100, text.scale = 2.5, nintersects = 40, order.by = "freq")
    dev.off()

    ### To find the interseected genes
    list_filter <- list(
        "BC_BCC_high" = get(filename2[1]),
        "BC_CM_high" = get(filename2[2]),
        "BC_CRC_high" = get(filename2[3]),
        "BC_EA_high" = get(filename2[4]),
        "BC_HCC_high" = get(filename2[5]),
        "BC_ICC_high" = get(filename2[6]),
        "BC_NSCLC_high" = get(filename2[7]),
        "BC_PDAC_high" = get(filename2[8]),
        "BC_RCC_high" = get(filename2[9]),
        "BC_SCC_high" = get(filename2[10]),
        "BC_UM_high" = get(filename2[11])
    )

    df2 <- data.frame(gene = unique(unlist(list_filter)))

    library(dplyr)
    df1 <- lapply(list_filter, function(x) {
        data.frame(gene = x)
    }) %>%
        bind_rows(.id = "path")

    df_int <- lapply(df2$gene, function(x) {
        # pull the name of the intersections
        intersection <- df1 %>%
            dplyr::filter(gene == x) %>%
            arrange(path) %>%
            pull("path") %>%
            paste0(collapse = "|")

        # build the dataframe
        data.frame(gene = x, int = intersection)
    }) %>%
        bind_rows()

    library(dplyr)
    final <- df_int %>%
        group_by(int) %>%
        dplyr::summarise(n = n()) %>%
        arrange(desc(n))

    df_int$int2 <- gsub("\\|", "__", df_int$int)

    write.table(df_int, paste(filepath2, "/Table/", cancername, "_df_int.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(final, paste(filepath2, "/Table/", cancername, "_final.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")

    ### We will be considering all the genes that are coming in 14,15 and 16 samples
    library(stringr)
    df_int$sample_num <- str_count(df_int$int, "_high")
    morethan_8 <- df_int[(df_int$sample_num >= 9), 1]
    morethan_7 <- df_int[(df_int$sample_num >= 8), 1]
    morethan_6 <- df_int[(df_int$sample_num >= 7), 1]
    morethan_5 <- df_int[(df_int$sample_num >= 6), 1]
    morethan_4 <- df_int[(df_int$sample_num >= 5), 1]
    morethan_3 <- df_int[(df_int$sample_num >= 4), 1]

    write.table(morethan_8, paste(filepath2, "/Table/", cancername, "_morethan_8.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(morethan_7, paste(filepath2, "/Table/", cancername, "_morethan_7.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(morethan_6, paste(filepath2, "/Table/", cancername, "_morethan_6.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(morethan_5, paste(filepath2, "/Table/", cancername, "_morethan_5.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(morethan_4, paste(filepath2, "/Table/", cancername, "_morethan_4.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(morethan_3, paste(filepath2, "/Table/", cancername, "_morethan_3.txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
    return(df_int)
}
