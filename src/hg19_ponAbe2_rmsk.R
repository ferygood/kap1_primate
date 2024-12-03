# this script prepare rmsk for hg19 and ponAbe2, because ponAbe2 rmsk is
# not in UCSC, so prepareRMSK() fail.

hg19 <- read.table("~/Downloads/hg19.fa.out", header = TRUE,
                   skip=1, fill=TRUE, stringsAsFactors = FALSE, row.names = NULL)

panTro4 <- read.table("~/Downloads/panTro4.fa.out", header=TRUE,
                      skip=1, fill=TRUE, stringsAsFactors = FALSE, row.names=NULL)

ponAbe2 <- read.table("~/Downloads/ponAbe3.fa.out", header = TRUE,
                      skip=1, fill=TRUE, stringsAsFactors = FALSE, row.names = NULL)

hg19_select <- hg19 %>%
    dplyr::select(c(sequence, begin, repeat., class.family))
colnames(hg19_select)[c(1, 2, 3, 4)] <- c("start", "end", "repName", "repClass")

panTro4_select <- panTro4 %>%
    dplyr::select(c(sequence, begin, repeat., class.family))
colnames(panTro4_select)[c(1, 2, 3, 4)] <- c("start", "end", "repName", "repClass")

ponAbe2_select <- ponAbe2 %>%
    dplyr::select(c(sequence, begin, repeat., class.family))
colnames(ponAbe2_select)[c(1, 2, 3, 4)] <- c("start", "end", "repName", "repClass")

hg19_len <- hg19_select %>%
    mutate(Len = end-start + 1) %>%
    group_by(repName, repClass) %>%
    summarise(rLen = mean(Len))

panTro4_len <- panTro4_select %>%
    mutate(Len = end-start + 1) %>%
    group_by(repName, repClass) %>%
    summarise(rLen = mean(Len))

ponAbe2_len <- ponAbe2_select %>%
    mutate(Len = end-start + 1) %>%
    group_by(repName, repClass) %>%
    summarise(cLen = mean(Len))

hm_oran_rmsk <- hg19_len %>%
    inner_join(ponAbe2_len[,c(1,3)], join_by(repName==repName))

chimp_oran_rmsk <- panTro4_len %>%
    inner_join(panTro4_len[,c(1,3)], join_by(repName==repName))
