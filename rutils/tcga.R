GeneSplitRegulon <- function(Genelist, Sep) {
    RegSplitted <- as.matrix(unlist(strsplit(as.character(Genelist), Sep)))

    return(RegSplitted)
}


TCGAanalyze_EA <- function(GeneName, RegulonList, TableEnrichment, EAGenes, GOtype,
                           FDRThresh = 0.01) {
    topPathways <- nrow(TableEnrichment)
    topPathways_tab <- matrix(0, 1, topPathways)
    topPathways_tab <- as.matrix(topPathways_tab)
    rownames(topPathways_tab) <- GeneName
    rownames(EAGenes) <- toupper(rownames(EAGenes))
    EAGenes <- EAGenes[!duplicated(EAGenes[, "ID"]), ]
    rownames(EAGenes) <- EAGenes[, "ID"]
    allgene <- EAGenes[, "ID"]
    current_pathway_from_EA <- as.matrix(TableEnrichment[, GOtype])
    TableNames <- gsub("David", "", paste("Top ", GOtype, " n. ",
        1:topPathways, " of ", topPathways,
        sep = ""
    ))
    colnames(topPathways_tab) <- TableNames
    topPathways_tab <- as.data.frame(topPathways_tab)
    table_pathway_enriched <- matrix(
        1, nrow(current_pathway_from_EA),
        8
    )
    colnames(table_pathway_enriched) <- c(
        "Pathway", "GenesInPathway",
        "Pvalue", "FDR", "CommonGenesPathway", "PercentPathway",
        "PercentRegulon", "CommonGeneSymbols"
    )
    table_pathway_enriched <- as.data.frame(table_pathway_enriched)
    for (i in 1:nrow(current_pathway_from_EA)) {
        table_pathway_enriched[i, "Pathway"] <- as.character(current_pathway_from_EA[i, ])
        if (nrow(TableEnrichment) == 589) {
            genes_from_current_pathway_from_EA <- GeneSplitRegulon(TableEnrichment[TableEnrichment[GOtype] ==
                as.character(current_pathway_from_EA[i, ]), ][
                ,
                "Molecules"
            ], ",")
        }
        else {
            genes_from_current_pathway_from_EA <- GeneSplitRegulon(TableEnrichment[TableEnrichment[GOtype] ==
                as.character(current_pathway_from_EA[i, ]), ][
                ,
                "Molecules"
            ], ", ")
        }
        genes_common_pathway_TFregulon <- as.matrix(intersect(
            toupper(RegulonList),
            toupper(genes_from_current_pathway_from_EA)
        ))
        if (length(genes_common_pathway_TFregulon) != 0) {
            current_pathway_commongenes_num <- length(genes_common_pathway_TFregulon)
            seta <- allgene %in% RegulonList
            setb <- allgene %in% genes_from_current_pathway_from_EA
            ft <- fisher.test(seta, setb)
            FisherpvalueTF <- ft$p.value
            table_pathway_enriched[i, "Pvalue"] <- as.numeric(FisherpvalueTF)
            if (FisherpvalueTF < 0.01) {
                current_pathway_commongenes_percent <- paste(
                    "(",
                    format((current_pathway_commongenes_num / length(genes_from_current_pathway_from_EA)) *
                        100, digits = 2), "%)"
                )
                current_pathway_commongenes_num_with_percent <- gsub(
                    " ",
                    "", paste(
                        current_pathway_commongenes_num,
                        current_pathway_commongenes_percent, "pv=",
                        format(FisherpvalueTF, digits = 2)
                    )
                )
                table_pathway_enriched[i, "CommonGenesPathway"] <- length(genes_common_pathway_TFregulon)
                table_pathway_enriched[i, "CommonGeneSymbols"] <- paste0(genes_common_pathway_TFregulon, collapse = ";")
                table_pathway_enriched[i, "GenesInPathway"] <- length(genes_from_current_pathway_from_EA)
                table_pathway_enriched[i, "PercentPathway"] <- as.numeric(table_pathway_enriched[
                    i,
                    "CommonGenesPathway"
                ]) / as.numeric(table_pathway_enriched[
                    i,
                    "GenesInPathway"
                ]) * 100
                table_pathway_enriched[i, "PercentRegulon"] <- as.numeric(table_pathway_enriched[
                    i,
                    "CommonGenesPathway"
                ]) / length(RegulonList) *
                    100
            }
        }
    }
    table_pathway_enriched <- table_pathway_enriched[order(table_pathway_enriched[
        ,
        "Pvalue"
    ], decreasing = FALSE), ]
    table_pathway_enriched <- table_pathway_enriched[table_pathway_enriched[
        ,
        "Pvalue"
    ] < 0.01, ]
    table_pathway_enriched[, "FDR"] <- p.adjust(table_pathway_enriched[
        ,
        "Pvalue"
    ], method = "fdr")
    table_pathway_enriched <- table_pathway_enriched[table_pathway_enriched[
        ,
        "FDR"
    ] < FDRThresh, ]
    table_pathway_enriched <- table_pathway_enriched[order(table_pathway_enriched[
        ,
        "FDR"
    ], decreasing = FALSE), ]

    table_pathway_enriched_filt <- table_pathway_enriched[table_pathway_enriched$FDR < FDRThresh, ]

    if (nrow(table_pathway_enriched) > 0) {
        tmp <- table_pathway_enriched
        tmp <- paste(tmp[, "Pathway"], "; FDR= ", format(tmp[
            ,
            "FDR"
        ], digits = 3), "; (ng=", round(tmp[, "GenesInPathway"]),
        "); (ncommon=", format(tmp[, "CommonGenesPathway"],
            digits = 2
        ), ")",
        sep = ""
        )
        tmp <- as.matrix(tmp)
        topPathways_tab <- topPathways_tab[, 1:nrow(table_pathway_enriched),
            drop = FALSE
        ]
        topPathways_tab[1, ] <- tmp
    }
    else {
        topPathways_tab <- NA
    }

    return(table_pathway_enriched_filt)
}

TCGAanalyze_EAcomplete <- function(TFname, RegulonList) {

    # This is a verification of the input
    # in case the List is like Gene|ID
    # we will get only the Gene
    local_split <- function(x) {
      unlist(purrr::map(x,
          .f = function(x) {
              stringr::str_split(x, ";", simplify = TRUE)
          }
      ))}
    if (all(grepl("\\|", RegulonList))) {
        RegulonList <- strsplit(RegulonList, "\\|")
        RegulonList <- unlist(lapply(RegulonList, function(x) x[1]))
    }

    print(paste(
        "I need about ", "1 minute to finish complete ",
        "Enrichment analysis GO[BP,MF,CC] and Pathways... "
    ))
    load(here::here("rutils", "TCGAbiolinks_sysdata.rda"))

    ResBP <- TCGAanalyze_EA(TFname, RegulonList, DAVID_BP_matrix,
        EAGenes,
        GOtype = "DavidBP"
    )
    ResBP <- local_split(ResBP$CommonGeneSymbols)

    print("GO Enrichment Analysis BP completed....done")
    ResMF <- TCGAanalyze_EA(TFname, RegulonList, DAVID_MF_matrix,
        EAGenes,
        GOtype = "DavidMF"
    )
    ResMF <- local_split(ResMF$CommonGeneSymbols)
    print("GO Enrichment Analysis MF completed....done")
    ResCC <- TCGAanalyze_EA(TFname, RegulonList, DAVID_CC_matrix,
        EAGenes,
        GOtype = "DavidCC"
    )
    ResCC <- local_split(ResCC$CommonGeneSymbols)
    print("GO Enrichment Analysis CC completed....done")
    ResPat <- TCGAanalyze_EA(TFname, RegulonList, listEA_pathways,
        EAGenes,
        GOtype = "Pathway"
    )
    ResPat <- local_split(ResPat$CommonGeneSymbols)
    print("Pathway Enrichment Analysis completed....done")

    ans <- list(ResBP = ResBP, ResMF = ResMF, ResCC = ResCC, ResPat = ResPat)
    return(ans)
}
