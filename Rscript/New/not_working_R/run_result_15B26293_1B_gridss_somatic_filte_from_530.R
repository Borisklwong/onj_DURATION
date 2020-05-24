> #!/usr/bin/env Rscript
  > #For 18_TGEG_530
  > library(argparser)
> argp = arg_parser("Filters a raw GRIDSS VCF into somatic call subsets.")
> argp = add_argument(argp, "--pondir", default=NA, help="Directory containing Panel Of Normal bed/bedpe used to filter FP somatic events. USer create_gridss_pon.R to generate the PON.")
> #BiocManager::install("BSgenome")
  > #BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
  > argp = add_argument(argp, "--ref", default="BSgenome.Hsapiens.UCSC.hg19", help="Reference genome to use. Must be a valid installed BSgenome package")
  > argp = add_argument(argp, "--input", help="GRIDSS VCF")
  > argp = add_argument(argp, "--output", help="High confidence somatic subset")
  > argp = add_argument(argp, "--pairoutput", help="High confidence somatic subset, paired reult")
  > argp = add_argument(argp, "--pairoutput_details", help="High confidence somatic subset, paired reult with details")
  > argp = add_argument(argp, "--fulloutput", help="Full call set excluding obviously germline call.")
  > argp = add_argument(argp, "--scriptdir", default=ifelse(sys.nframe() == 0, "./", dirname(sys.frame(1)$ofile)), help="Path to libgridss_DURATION.R script")
  > argp = add_argument(argp, "--gc", flag=TRUE, help="Perform garbage collection after freeing of large objects. ")
  > argv = parse_args(argp, argv=c(
    +   "--input", "/stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_sort_mark.sv.vcf",
    +   "--output", "/stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_postfilter.sv.vcf",
    +   "--pairoutput", "/stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_postfilter_paired.sv.bedpe",
    +   "--pairoutput_details", "/stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_postfilter_paired.csv",
    +   "-f", "/stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_postfilter_full.sv.vcf",
    +   "-p", "/stornext/Home/data/allstaff/w/wong.b/R/GRIDSS_PON_3792v1",
    +   "--scriptdir", "/stornext/Home/data/allstaff/w/wong.b/R/R_script",
    +   "--gc"))
  > 
    > 
    > if (!file.exists(argv$input)) {
      +   msg = paste(argv$input, "not found")
      +   write(msg, stderr())
      +   print(argp)
      +   stop(msg)
      + }
  > if (is.na(argv$pondir)) {
    +   argv$pondir = NULL
    + } else if (!dir.exists(argv$pondir)) {
      +   msg = paste(argv$pondir, "not found")
      +   write(msg, stderr())
      +   print(argp)
      +   stop(msg)
      + }
  > refgenome=eval(parse(text=paste0("library(", argv$ref, ")\n", argv$ref)))
  > 
    > library(tidyverse)
  > library(readr)
  > library(stringr)
  > libgridssfile = paste0(argv$scriptdir, "/", "libgridss_DURATION.R")
  > if (file.exists(libgridssfile)) {
    +   tmpwd = getwd()
    +   setwd(argv$scriptdir)
    +   source("libgridss_DURATION.R")
    +   setwd(tmpwd)
    + } else {
      +   msg = paste("Could not find libgridss_DURATION.R in", argv$scriptdir, " - please specify a --scriptdir path to a directory containing the required scripts")
      +   write(msg, stderr())
      +   print(argp)
      +   stop(msg)
      + }
  > 
    > 
    > # Filter to somatic calls
    > write(paste(Sys.time(), "Reading", argv$input), stderr())
  2020-05-24 21:16:38 Reading /stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_sort_mark.sv.vcf
  > raw_vcf = readVcf(argv$input)
  > tumourordinal = 1
  > # hard filter variants that are too low quality
    > full_vcf = raw_vcf[VariantAnnotation::fixed(raw_vcf)$QUAL > 350 &
                           +                      !str_detect(as.character(seqnames(raw_vcf)), "Un") &
                           +                      !str_detect(as.character(seqnames(raw_vcf)), "gl")]
  > rm(raw_vcf)
  > if (argv$gc) { gc() }
  used  (Mb) gc trigger   (Mb)  max used   (Mb)
  Ncells  9978937 533.0   28403569 1517.0  44518326 2377.6
  Vcells 48530884 370.3  179894727 1372.5 554950879 4234.0
  > # hard filter unpaired breakpoints (caused by inconsistent scoring across the two breakends)
    > full_vcf = full_vcf[is.na(info(full_vcf)$PARID) | info(full_vcf)$PARID %in% names(full_vcf)]
  > #full_vcf = align_breakpoints(full_vcf)
    > # Add header fields
    > full_vcf = addVCFHeaders(full_vcf)
  > 
    > info(full_vcf)$BPI_AF = rep("", length(full_vcf))
  > filters = rep("BE", length(full_vcf))
  > names(filters) = names(full_vcf)
  > 
    > 
    > write(paste(Sys.time(), "Parsing breakpoints", argv$input), stderr())
  2020-05-24 21:17:20 Parsing breakpoints /stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_sort_mark.sv.vcf
  > bpgr = breakpointRanges(full_vcf, unpartneredBreakends=FALSE)
  > bpgr$vcfId = names(bpgr)
  > write(paste(Sys.time(), "Calculating breakpoint VAF", argv$input), stderr())
  2020-05-24 21:17:21 Calculating breakpoint VAF /stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_sort_mark.sv.vcf
  > bpgr$af = round(gridss_bp_af(bpgr, full_vcf, tumourordinal), 5)
  > bpgr$af_str = paste(bpgr$af, partner(bpgr)$af, sep=",")
  > if (length(bpgr) > 0) {
    +   info(full_vcf[names(bpgr)])$BPI_AF = bpgr$af_str
    + }
  > write(paste(Sys.time(), "Filtering breakpoints", argv$input), stderr())
  2020-05-24 21:17:21 Filtering breakpoints /stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_sort_mark.sv.vcf
  > bpfiltered = gridss_breakpoint_filter(bpgr, full_vcf, bsgenome=refgenome, pon_dir=argv$pondir)
  2020-05-24 21:17:21 Loading from cache  /stornext/Home/data/allstaff/w/wong.b/R/GRIDSS_PON_3792v1/gridss_pon_breakpoint.bedpe.cached.parsed.rds
  > filters[names(bpgr)] = bpfiltered
  > if (argv$gc) { gc() }
  used  (Mb) gc trigger   (Mb)  max used   (Mb)
  Ncells  9490194 506.9   36396129 1943.8  44518326 2377.6
  Vcells 34920536 266.5  358584778 2735.8 554950879 4234.0
  > 
    > # shadow breakpoint removed due to initial mapq20 filter reducing FP rate
    > # bpfiltered = .addFilter(bpfiltered, "shadow", is_shadow_breakpoint(bpgr, begr, full_vcf))
    > 
    > #bpfiltered = .addFilter(bpfiltered, "LOW_Qual", bpgr$QUAL < gridss.min_qual)
    > #som_llr = gridss_breakpoint_somatic_llr(full_vcf, bpgr=bpgr, contamination_rate=gridss.allowable_normal_contamination)
    > 
    > # - filter to only decent length assemblies?
    > #begr$calls_1k_window = countOverlaps(begr, rowRanges(full_vcf), ignore.strand=TRUE, maxgap=1000)
    > 
    > # Remove very hard filtered variants
    > full_vcf = full_vcf[passes_very_hard_filters(filters)]
  > unpaired_breakpoint = !is.na(info(full_vcf)$PARID) & !(info(full_vcf)$PARID %in% names(full_vcf))
  > full_vcf = full_vcf[!unpaired_breakpoint]
  > filters = filters[passes_very_hard_filters(filters)]
  > filters = filters[!unpaired_breakpoint]
  > rm(unpaired_breakpoint)
  > 
    > vcf = full_vcf[passes_soft_filters(filters)]
  > vcf = vcf[is.na(info(vcf)$PARID) | info(vcf)$PARID %in% names(vcf)]
  > bpgr = bpgr[names(bpgr) %in% names(vcf)]
  > if (argv$gc) { gc() }
  used  (Mb) gc trigger   (Mb)  max used   (Mb)
  Ncells  9490467 506.9   29116904 1555.1  44518326 2377.6
  Vcells 34920446 266.5  286867823 2188.7 554950879 4234.0
  > 
    > 
    > 
    > 
    > 
    > 
    > 
    > 
    >  
    > 
    > #####################
  > #### here on is optional
    > 
    > 
    > # write(paste(Sys.time(),"Calculating transitive links", argv$input), stderr())
    > # # transitive calling
    > # transitive_df = transitive_calls(vcf, bpgr, report="max2") %>%
    > #   # only make transitive calls were we actually know the path
    > #   filter(!has_multiple_paths) %>%
    > #   mutate(type="transitive")
    > # # now we filter imprecise variants,# can be skip
    > # is_imprecise = !(is.na(info(vcf)$IMPRECISE) | !info(vcf)$IMPRECISE) |
    > #   !((!is.na(info(vcf)$PARID) & info(vcf)$ASSR + info(vcf)$SR + info(vcf)$IC > 0) |
    > #       (is.na(info(vcf)$PARID) & info(vcf)$BASSR + info(vcf)$BSC > 0))
    > # filters[names(vcf)[is_imprecise]] = paste0(filters[names(vcf)[is_imprecise]], ";imprecise")
    > # filters[transitive_df$linked_by] = paste0(filters[transitive_df$linked_by], ";transitive")
    > # ############################################################
  > # vcf = full_vcf[passes_soft_filters(filters)]
    > # vcf = vcf[is.na(info(vcf)$PARID) | info(vcf)$PARID %in% names(vcf)]
    > # bpgr = bpgr[names(bpgr) %in% names(vcf)]
    > # begr = begr[names(begr) %in% names(vcf)]
    > # if (argv$gc) { gc() }
    > 
    > asm_linked_df = NULL
  > write(paste(Sys.time(),"Calculating assembly links", argv$input), stderr())
  2020-05-24 21:19:03 Calculating assembly links /stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_sort_mark.sv.vcf
  > # Assembly-based event linking
    > if (length(vcf) > 0) {
      +   asm_linked_df = linked_assemblies(vcf) %>%
        +     mutate(type="asm")
      + }
  > 
    > link_df = bind_rows(asm_linked_df) %>% #, transitive_df) %>%
    +   mutate(linking_group=str_replace(linked_by, "/.*$", "")) %>%
    +   mutate(pass=passes_final_filters(vcf[vcfId])) %>%
    +   group_by(linking_group) %>%
    +   mutate(pass=any(pass)) %>%
    +   ungroup() %>%
    +   filter(pass)
  Error in stri_replace_first_regex(string, pattern, fix_replacement(replacement),  : 
                                      object 'linked_by' not found
                                    > 
                                      > 
                                      > 
                                      > # Inversion linkage
                                      > write(paste(Sys.time(),"Calculating simple inversions", argv$input), stderr())
                                    2020-05-24 21:19:03 Calculating simple inversions /stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_sort_mark.sv.vcf
                                    > inv_link_df = linked_by_simple_inversion_classification(bpgr) %>%
                                      +   group_by(linked_by) %>%
                                      +   mutate(pass=passes_final_filters(vcf[vcfId])) %>%
                                      +   mutate(pass=any(pass)) %>%
                                      +   ungroup() %>%
                                      +   filter(pass) %>%
                                      +   mutate(type="inv")
                                    > # Deletion bridge linkage
                                      > # TODO: do we want to do this?
                                      > # I'm suspicious of the model used in ChainFinder PMC3690918
                                      > # Notably: I'm suspicous that "repair with major DNA loss" is actually a thing
                                      > # given the catastrophic nature of chromo*, a more reasonable explaination is
                                      > # an additional DSB with the subsequent loss of that DNA fragment.
                                      > # Given the focal nature of chromoplexy, ChainFinder works because it just
                                      > # finds the focal events, not because the model is correct.
                                      > # TODO: show this by modelling additional focal DSBs
                                      > write(paste(Sys.time(),"Calculating dsb links", argv$input), stderr())
                                    2020-05-24 21:19:03 Calculating dsb links /stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_sort_mark.sv.vcf
                                    > dsb_link_df = linked_by_dsb(bpgr) %>%
                                      +   group_by(linked_by) %>%
                                      +   mutate(pass=passes_final_filters(vcf[vcfId])) %>%
                                      +   mutate(pass=any(pass)) %>%
                                      +   ungroup() %>%
                                      +   filter(pass) %>%
                                      +   mutate(type="dsb")
                                    Warning messages:
                                      1: In min(distance) : no non-missing arguments to min; returning Inf
                                    2: In min(qqual) : no non-missing arguments to min; returning Inf
                                    3: In min(distance) : no non-missing arguments to min; returning Inf
                                    4: In min(squal) : no non-missing arguments to min; returning Inf
                                    > 
                                      > write(paste(Sys.time(),"Removing duplicated/conflicting links", argv$input), stderr())
                                    2020-05-24 21:19:04 Removing duplicated/conflicting links /stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_sort_mark.sv.vcf
                                    > # linking priorities:
                                      > # - asm independent of other linkages
                                      > # - transitive independent of other linkages
                                      > # - ins, inv, dsb linkages
                                      > event_link_df = bind_rows(
                                        +   inv_link_df,
                                        +   dsb_link_df) %>%
                                      +   dplyr::select(vcfId, linked_by) %>%
                                      +   mutate(
                                        +     QUAL=rowRanges(vcf)[vcfId]$QUAL,
                                        +     hasPolyA=str_detect(rowRanges(vcf[vcfId])$ALT, "A{16}")) %>%
                                      +   group_by(linked_by) %>%
                                      +   # filter events where supporting fragment counts differ by too much
                                      +   mutate(
                                        +     max_supporting_fragment_count = max(ifelse(is.na(info(full_vcf[vcfId])$PARID), info(full_vcf[vcfId])$BVF, info(full_vcf[vcfId])$VF)),
                                        +     min_supporting_fragment_count = min(ifelse(is.na(info(full_vcf[vcfId])$PARID), info(full_vcf[vcfId])$BVF, info(full_vcf[vcfId])$VF)),
                                        +     hasPolyA=any(hasPolyA)
                                        +     ) %>%
                                      +   filter(min_supporting_fragment_count >= gridss.min_rescue_portion * max_supporting_fragment_count | hasPolyA)
                                    Warning messages:
                                      1: In max(ifelse(is.na(info(full_vcf[vcfId])$PARID), info(full_vcf[vcfId])$BVF,  :
                                                         no non-missing arguments to max; returning -Inf
                                                       2: In min(ifelse(is.na(info(full_vcf[vcfId])$PARID), info(full_vcf[vcfId])$BVF,  :
                                                                          no non-missing arguments to min; returning Inf
                                                                        > 
                                                                          > write(paste(Sys.time(),"Calculating final linkage annotation", argv$input), stderr())
                                                                        2020-05-24 21:19:04 Calculating final linkage annotation /stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_sort_mark.sv.vcf
                                                                        > # Only keep the best QUAL event linkage
                                                                          > event_link_df = event_link_df %>%
                                                                          +   group_by(linked_by) %>%
                                                                          +   mutate(linkQUAL = pmin(QUAL)) %>%
                                                                          +   group_by(vcfId) %>%
                                                                          +   filter(QUAL == linkQUAL) %>%
                                                                          +   group_by(linked_by)
                                                                        > # Don't event link to PON filtered variants
                                                                          > event_link_df = event_link_df %>%
                                                                          +   filter(!str_detect(filters[vcfId], "PON"))
                                                                        > # Fix up pairing
                                                                          > event_link_df = event_link_df %>%
                                                                          +   filter(n() == 2) %>%
                                                                          +   ungroup()
                                                                        > 
                                                                          > # include both breakends of any linked breakpoints
                                                                          > # as linkage can be breakend specific (e.g. assembly, bpbeins)
                                                                          > link_rescue = bind_rows(link_df, event_link_df) %>% pull(vcfId) %>% unique()
                                                                        Error in dots_values(...) : object 'link_df' not found
                                                                        > link_rescue = c(link_rescue, bpgr[link_rescue[link_rescue %in% names(bpgr)]]$partner)
                                                                        Error: object 'link_rescue' not found
                                                                        > 
                                                                          > # Note that we don't rescue equivalent events
                                                                          > eqv_link_df = linked_by_equivalent_variants(full_vcf, as(rbind(as.data.frame(bpgr)), "GRanges"), bsgenome=refgenome) %>%
                                                                          +   filter(passes_final_filters(vcf[vcfId]) | vcfId %in% link_rescue) %>%
                                                                          +   group_by(linked_by) %>%
                                                                          +   filter(n() == 2) %>%
                                                                          +   ungroup() %>%
                                                                          +   mutate(type="eqv")
                                                                        > 
                                                                          > link_summary_df = bind_rows(link_df, event_link_df, eqv_link_df) %>%
                                                                          +   group_by(vcfId) %>%
                                                                          +   summarise(linked_by=paste0(linked_by, collapse=","))
                                                                        Error in dots_values(...) : object 'link_df' not found
                                                                        > 
                                                                          > # Add linking information
                                                                          > info(full_vcf)$LOCAL_LINKED_BY = rep("", length(full_vcf))
                                                                        > info(full_vcf)$REMOTE_LINKED_BY = rep("", length(full_vcf))
                                                                        > info(full_vcf[link_summary_df$vcfId])$LOCAL_LINKED_BY = link_summary_df$linked_by
                                                                        Error: object 'link_summary_df' not found
                                                                        > info(full_vcf[!is.na(info(full_vcf)$PARID)])$REMOTE_LINKED_BY = info(full_vcf[info(full_vcf[!is.na(info(full_vcf)$PARID)])$PARID])$LOCAL_LINKED_BY
                                                                        > 
                                                                          > write(paste(Sys.time(),"Final qual filtering ", argv$output), stderr())
                                                                        2020-05-24 21:19:04 Final qual filtering  /stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_postfilter.sv.vcf
                                                                        > # final qual filtering
                                                                          > fails_qual_without_rescue = !passes_final_filters(full_vcf, include.existing.filters=FALSE) & !(names(full_vcf) %in% link_rescue)
                                                                        Error in names(full_vcf) %in% link_rescue : 
                                                                          object 'link_rescue' not found
                                                                        > filters[names(full_vcf)[fails_qual_without_rescue]] = paste0(filters[names(full_vcf)[fails_qual_without_rescue]], ";qual")
                                                                        Error in paste0(filters[names(full_vcf)[fails_qual_without_rescue]], ";qual") : 
                                                                          object 'fails_qual_without_rescue' not found
                                                                        > 
                                                                          > 
                                                                          > 
                                                                          > 
                                                                          > 
                                                                          > 
                                                                          > ################
                                                                        > # Write outputs
                                                                          > VariantAnnotation::fixed(full_vcf)$FILTER = ifelse(str_remove(filters, "^;") == "", "PASS", str_remove(filters, "^;"))
                                                                        > if (argv$gc) { gc() }
                                                                        used  (Mb) gc trigger   (Mb)  max used   (Mb)
                                                                        Ncells  9488694 506.8   29116904 1555.1  44518326 2377.6
                                                                        Vcells 34920013 266.5  229494259 1751.0 554950879 4234.0
                                                                        > if (!is.na(argv$output)) {
                                                                          +   write(paste(Sys.time(),"Writing ", argv$output), stderr())
                                                                          +   vcf = full_vcf[passes_soft_filters(filters)]
                                                                          +   vcf = vcf[is.na(info(vcf)$PARID) | info(vcf)$PARID %in% names(vcf)]
                                                                          +   writeVcf(vcf, argv$output)
                                                                          + }
                                                                        2020-05-24 21:19:05 Writing  /stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_postfilter.sv.vcf
                                                                        > 
                                                                          > 
                                                                          > out_vcf.file = argv$output
                                                                        > out_bpgr = breakpointRanges(VariantAnnotation::readVcf(out_vcf.file))
                                                                        > #Only keep the interchromosomal SVs
                                                                          > out_bpgr_inter = out_bpgr[is.na(out_bpgr$svLen)]
                                                                        > #Pair the results
                                                                          > out_pairgr = breakpointgr2pairs(out_bpgr_inter)
                                                                        > #sort by score
                                                                          > out_pairgr_sort = out_pairgr[order(out_pairgr@first@elementMetadata@listData[["QUAL"]], decreasing = TRUE),]
                                                                        > write.csv(out_pairgr_sort, argv$pairoutput_details)
                                                                        > rtracklayer::export(out_pairgr_sort, con=argv$pairoutput)
                                                                        > 
                                                                          > 
                                                                          > 
                                                                          > 
                                                                          > 
                                                                          > if (!is.na(argv$fulloutput)) {
                                                                            +   write(paste(Sys.time(),"Writing ", argv$fulloutput), stderr())
                                                                            +   vcf = full_vcf[passes_very_hard_filters(filters)]
                                                                            +   vcf = vcf[is.na(info(vcf)$PARID) | info(vcf)$PARID %in% names(vcf)]
                                                                            +   writeVcf(vcf, argv$fulloutput)
                                                                            + }
                                                                        2020-05-24 21:19:06 Writing  /stornext/Home/data/allstaff/w/wong.b/R/R_script/MESO/15B26293_1B_script/15B26293_1B_postfilter_full.sv.vcf
                                                                        > 
                                                                          > 
                                                                          > ##Ploting circos plot
                                                                          > library(circlize)
                                                                        > #change the format of the result from Granges into bedpe
                                                                          > out_bpgr_inter_bedpe = breakpointgr2bedpe(out_bpgr_inter)
                                                                        > plot_title=paste(gsub(pattern = "_merged_mark.sv.vcf", "", basename(argv$input)),".jpeg",sep = "")
                                                                        > jpeg(file=plot_title)
                                                                        > circos.initializeWithIdeogram()
                                                                        > title(gsub(pattern = "_merged_mark.sv.vcf", "", basename(argv$input)))
                                                                        > circos.genomicLink(out_bpgr_inter_bedpe[,1:3], out_bpgr_inter_bedpe[,4:6])
                                                                        > dev.off()
                                                                        null device 
                                                                        1 
                                                                        > git add -A
                                                                        Error: unexpected symbol in "git add"