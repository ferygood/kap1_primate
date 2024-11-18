#' Estimate Repeated Measures Correlation comparing orthologs and TEs
#' @description
#' This function is a modified version of TEKRABber::corrOrthologTE(). It is
#' designed to handle repeat measurements from the same biological sample. For
#' example, you sample three times from an A cell, and named your samples
#' as A_1, A_2, A_3. This function leverage rmcorr package. You can find more
#' details in its paper (DOI:10.3389/fpsyg.2017.00456)
#'
#' @param geneInput
#' @param teInput
#' @param padjMethod
#' @param numCore
#' @param sample_meta
#' @param fileDir
#' @param fileName
#'
#' @return
#' @export
#'
#' @examples
corrOrthologTE_rmcorr <- function(geneInput, teInput, padjMethod="fdr",
                                  numCore=1, sample_meta, fileDir=NULL,
                                  fileName="TEKRABber_geneTE_RmCorr.csv"){

    # register backend
    registerDoParallel(numCore)

    # create empty dataframe
    df.corr <- data.frame(matrix(ncol=4, nrow=0))
    colnames(df.corr) <- c("geneName", "teName", "coef", "pvalue")

    # calculate correlation with parallel
    df.corr <- foreach(i=1:nrow(geneInput), .combine=rbind) %dopar% {
        foreach(j=1:nrow(teInput), .combine=rbind) %dopar% {
            num1 <- as.numeric(geneInput[i, ])
            num2 <- as.numeric(teInput[j, ])

            # create input for rmc
            rmc_input <- data.frame(
                sample_ID = sample_meta,
                gene = num1,
                te = num2
            )

            # calculate repeat correlation
            rmc <- rmcorr::rmcorr(
                participant = sample_ID,
                measure1 = gene,
                measure2 = te,
                dataset = rmc_input
            )

            rbind(
                data.frame(geneName=rownames(geneInput)[i],
                           teName=rownames(teInput)[j],
                           coef=rmc$r,
                           pvalue=rmc$p),
                df.corr
            )
        }
    }

    stopImplicitCluster()

    # calculate p-adjusted value
    df.corr$padj <- p.adjust(
        df.corr$pvalue,
        method=padjMethod
    )

    rownames(df.corr) <- seq_len(nrow(df.corr))

    if (!is.null(fileDir)){
        dir.create(fileDir)
        write.table(df.corr, file = file.path(fileDir, fileName), sep=",")
    }

    df.corr

}
