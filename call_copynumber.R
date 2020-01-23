# This pipeline is the default QDNAseq, which does not use the matched normal

library(QDNAseq)
library(Biobase)
library(png)
library(grid)
library(gridExtra)
library(gtools)
library(ggplot2)
args = commandArgs(T)
samplename = args[1]
tumourbam = args[2]
normalname = args[3]
normalbam = args[4]

if (length(args) > 4) {
	sex = args[5]
} else {
	sex = NULL
}
sex = "female"

correctReplication = F
segmentation_gamma = 30
segmentation_kmin = 3

# parameters for testing for copy number alterations
logr_minimum_deviation = 0.05 # minimum logr has to deviate from 0 to be considered for significance testing
max_num_bins_test = 500 # maximum number of bins a segment can contain before downsampling is done - for large numbers of bins a very small deviation will always be significant
test_significance_threshold = 0.05 # significance threshold

load("precalculated_windows/QNDAseq_bins50.RData")
readCounts_tumour <- binReadCounts(bins, bamfiles=tumourbam)
readCounts_normal <- binReadCounts(bins, bamfiles=normalbam)

# determine the sex
#if (is.null(sex)) {
#	y_counts = assayDataElement(readCounts_normal, "counts")[bins$chromosome=="Y"]
#	if (median(y_counts) > 1) {
#		sex = "male"
#	} else {
		sex = "female"
		# remove Y chromosome
		readCounts_tumour = new('QDNAseqReadCounts', bins=bins[bins$chromosome!="Y",,drop=F], counts=assayDataElement(readCounts_tumour, "counts")[bins$chromosome!="Y",,drop=F], phenodata=phenoData(readCounts_tumour))
		readCounts_normal = new('QDNAseqReadCounts', bins=bins[bins$chromosome!="Y",,drop=F], counts=assayDataElement(readCounts_normal, "counts")[bins$chromosome!="Y",,drop=F], phenodata=phenoData(readCounts_normal))
#	}
#}


normaliseReadCounts = function(readCounts_tumour, readCounts_normal, sample_index=1, minReadsThreshold=10, correctReplication=F) {
	tumour_counts = assayDataElement(readCounts_tumour, "counts")[,sample_index]
	normal_counts = assayDataElement(readCounts_normal, "counts")[,sample_index]

	normal_counts[normal_counts < minReadsThreshold] = NA
	
	tumour_r = tumour_counts/normal_counts
	tumour_r[!is.finite(tumour_r)] = NA
	tumour_logr = log2(tumour_r / mean(tumour_r, na.rm=T))
	tumour_logr[!is.finite(tumour_logr)] = NA
	tumour_logr = matrix(tumour_logr, ncol=1)
	dimnames(tumour_logr) = dimnames(assayDataElement(readCounts_tumour, "counts"))
	
	
	# create a QDNAseq object
	res = new("QDNAseqCopyNumbers", bins=featureData(readCounts_tumour), copynumber=tumour_logr,phenodata=phenoData(readCounts_tumour))
	
	# correct for gc and mappability
	counts = tumour_logr
	gc <- round(fData(res)$gc)
	mappability <- round(fData(res)$mappability)

	if (correctReplication) {
		replication <- round(fData(res)$replication)
		condition <- QDNAseq:::binsToUse(res) & !is.na(gc) & !is.na(mappability) & is.finite(replication)
		smoothT <- loess(counts[,sample_index]~gc*mappability*replication)
	} else {	
		condition <- QDNAseq:::binsToUse(res) & !is.na(gc) & !is.na(mappability)
		smoothT <- loess(counts[,sample_index]~gc*mappability)
	}
	fit = matrix(NA, ncol=1, nrow=nrow(res))
	dimnames(fit) = dimnames(counts)
	fit[names(smoothT$residuals),sample_index] = 2^smoothT$residuals
	assayDataElement(res, "copynumber") = fit
	QDNAseq:::binsToUse(res) <- QDNAseq:::binsToUse(res) & !is.na(rowMeans(fit))
	return(res)
}

segmentBinsPCF = function(copyNumbersSmooth, segmentation_gamma, segmentation_kmin) {
  condition = QDNAseq:::binsToUse(copyNumbersSmooth)	
  logr = log2(assayDataElement(copyNumbersSmooth, "copynumber")[condition,1])
  
  chroms = unique(fData(copyNumbersSmooth)$chromosome)
  
  logr_segm = c()
  for (chrom in chroms) {
    logr_chrom = logr[(fData(copyNumbersSmooth)$chromosome==chrom)[condition]]
    sdev = Battenberg:::getMad(logr_chrom, k=25)
    logr_segm = c(logr_segm, Battenberg:::selectFastPcf(logr_chrom, segmentation_kmin, segmentation_gamma*sdev, T)$yhat)
  }
  
  temp = assayDataElement(copyNumbersSmooth, "copynumber")
  temp[,1] = NA
  temp[condition,1] = 2^logr_segm
  assayDataElement(copyNumbersSmooth, "segmented") = temp
  
  return(copyNumbersSmooth)
}

normalisedCounts = normaliseReadCounts(readCounts_tumour, readCounts_normal)
copyNumbersSmooth <- smoothOutlierBins(normalisedCounts)
# copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented = segmentBinsPCF(copyNumbersSmooth, segmentation_gamma=segmentation_gamma, segmentation_kmin=segmentation_kmin)
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
copyNumbersCalled <- callBins(copyNumbersSegmented)

# make result figure of adapted QDNAseq output
#png(paste0("output/", samplename, "/", samplename, "_figure.png"), width=1000, height=400)
#plot(copyNumbersCalled)
#dev.off()

###############################################################################################
# make raw data figures
###############################################################################################
processRawDataDefault = function(readCounts) {
	copyNumbers <- correctBins(readCounts, correctReplication=T, method="median")
	copyNumbersNormalized <- normalizeBins(copyNumbers)
	copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
	copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
	copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
	return(copyNumbersSegmented)
}

tumour_raw = processRawDataDefault(readCounts_tumour)
normal_raw = processRawDataDefault(readCounts_normal)

#png(paste0("output/", samplename, "/", samplename, "_tumour_raw.png"), width=1000, height=400)
#plot(tumour_raw)
#dev.off()

#png(paste0("output/", samplename, "/", samplename, "_normal_raw.png"), width=1000, height=400)
#plot(normal_raw)
#dev.off()

png(paste0("output/", samplename, "/", samplename, "_rawdata.png"), width=1000, height=800)
par(mfrow=c(2,1))
plot(tumour_raw)
plot(normal_raw)
dev.off()


###############################################################################################
# make posthoc plot and segmentation output
###############################################################################################

get_segments = function(x) {
  ns <- asNamespace("CGHbase")
  .makeSegments <- get(".makeSegments", envir=ns, mode="function")
  
  condition <- QDNAseq:::binsToUse(x)	
  all.chrom <- QDNAseq:::chromosomes(x)
  if (is.integer(all.chrom)) # when x is a cghRaw, cghSeg, or cghCall object
    all.chrom <- as.character(all.chrom)
  chrom <- all.chrom[condition]
  
  condition <- QDNAseq:::binsToUse(x)
  copynumber <- as.data.frame(Biobase::assayDataElement(x, "copynumber")[condition, , drop=FALSE])
  
  # assuming a single sample here
  segmented <- Biobase::assayDataElement(x, "segmented")[condition, 1]
  segmented <- QDNAseq:::log2adhoc(segmented)
  segment <- as.data.frame(.makeSegments(segmented, chrom))
  
  copynumber$segment_values = NA
  # segmented$p_value_gain = NA
  # segmented$p_value_loss = NA
  for (i in seq_len(nrow(segment))) {
    copynumber$segment_values[segment$start[i]:segment$end[i]] = segment$values[i]
    # logr_values = log2(copynumber[segment$start[i]:segment$end[i], 1])
    # 
    # num_snps_sample = 500
    # offset = 0.01
    # 
    # copynumber$p_value_loss[i] = t.test(sample(logr_values, num_snps_sample), rnorm(n=num_snps_sample, mean=0-offset, sd=sd(logr_values)), alternative="less")$p.value
    # copynumber$p_value_gain[i] = t.test(sample(logr_values, num_snps_sample), rnorm(n=num_snps_sample, mean=0+offset, sd=sd(logr_values)), alternative="greater")$p.value
  }
  
  # copynumber$p_value_gain_adj = p.adjust(copynumber$p_value_gain, method = "bonferroni")
  # copynumber$p_value_loss_adj = p.adjust(copynumber$p_value_loss, method = "bonferroni")
  
  
  colnames(copynumber)[1] = "values"
  copynumber$logr_values = QDNAseq:::log2adhoc(copynumber$values)
  copynumber$chrom = unlist(lapply(rownames(copynumber), function(x) unlist(stringr::str_split(x, ":"))[1]))
  copynumber$start = as.numeric(unlist(lapply(rownames(copynumber), function(x) unlist(stringr::str_split(unlist(stringr::str_split(x, ":"))[2], "-"))[1])))
  copynumber$end = as.numeric(unlist(lapply(rownames(copynumber), function(x) unlist(stringr::str_split(unlist(stringr::str_split(x, ":"))[2], "-"))[2])))
  
  copynumber$chrom = factor(copynumber$chrom, levels=gtools::mixedsort(unique(copynumber$chrom)))
  
  return(copynumber)
}
copynumber = get_segments(copyNumbersCalled)


get_segmentation_df = function(copynumber, max_num_bins_test, logr_minimum_deviation, test_significance_threshold) {
  segments = rle(copynumber$segment_values)
  segmentation = data.frame()
  for (i in seq_along(segments$lengths)) {
    if (i==1) {
      curr_segment_start_index = 1
    } else {
      segment_end_index = cumsum(segments$lengths[1:(i-1)])
      curr_segment_start_index = segment_end_index[length(segment_end_index)] + 1
    }
    
    if (i==length(segments$lengths)) {
      curr_segment_end_index = nrow(copynumber)
    } else {
      segment_end_index = cumsum(segments$lengths[1:i])
      curr_segment_end_index = segment_end_index[length(segment_end_index)]
    }
    
    # test for significance of gain and loss separately
    logr_values = log2(copynumber$values[curr_segment_start_index:curr_segment_end_index])
    num_snps_sample = ifelse(length(logr_values) > max_num_bins_test, max_num_bins_test, length(logr_values))
    p_value_gain = t.test(sample(logr_values, num_snps_sample), rnorm(n=num_snps_sample, mean=0+logr_minimum_deviation, sd=sd(logr_values)), alternative="greater")$p.value
    p_value_loss = t.test(sample(logr_values, num_snps_sample), rnorm(n=num_snps_sample, mean=0-logr_minimum_deviation, sd=sd(logr_values)), alternative="less")$p.value
    
    segmentation = rbind(segmentation, data.frame(chrom=copynumber$chrom[curr_segment_start_index],
                                                  start=copynumber$start[curr_segment_start_index], end=copynumber$end[curr_segment_end_index],
                                                  value=copynumber$segment_values[curr_segment_start_index], 
                                                  p_value_gain=p_value_gain,
                                                  p_value_loss=p_value_loss, stringsAsFactors=F))
  }
  
  # adjust for multiple testing
  segmentation$p_value_gain_adj = p.adjust(segmentation$p_value_gain, method = "bonferroni")
  segmentation$p_value_loss_adj = p.adjust(segmentation$p_value_loss, method = "bonferroni")
  # classify segments
  segmentation$classification = "normal"
  is_gain = segmentation$p_value_gain_adj < test_significance_threshold
  is_loss = segmentation$p_value_loss_adj < test_significance_threshold
  segmentation$classification[is_gain] = "gain"
  segmentation$classification[is_loss] = "loss"
  # check that no segments are both a gain and a loss - this should never occur
  if (any(is_gain & is_loss)) {
    segmentation$classification[is_gain & is_loss] = "error"
    print("Found segment(s) that is both gain and loss:")
    print(segmentation[is_gain & is_loss,])
  }
  segmentation$len = segmentation$end - segmentation$start
  return(segmentation)
}
segmentation = get_segmentation_df(copynumber, max_num_bins_test=max_num_bins_test, logr_minimum_deviation=logr_minimum_deviation, test_significance_threshold=test_significance_threshold)

get_purity_estimates = function(segmentation, sex) {
  segmentation$purity = NA
  # estimate purity based on calls from qdnaseq?
  chrom_lengths = data.frame()
  for (chrom in unique(copynumber$chrom)) {
    pos_start = min(copynumber[copynumber$chrom==chrom, "start"])
    pos_end = max(copynumber[copynumber$chrom==chrom, "end"])
    chrom_lengths = rbind(chrom_lengths, data.frame(chrom=chrom, start=pos_start, end=pos_end))
  }
  chrom_lengths$len = chrom_lengths$end - chrom_lengths$start
  chrom_lengths$frac = (chrom_lengths$len/1000) / sum(chrom_lengths$len / 1000)
  
  # obtain a ploidy, assuming this sample is normal
  if (sex=="male") {
    normal_ploidy = sum(c(chrom_lengths$len[1:19]/1000*2, chrom_lengths$len[20]/1000*1, chrom_lengths$len[21]/1000*1)) / sum(chrom_lengths$len / 1000)
  } else {
    normal_ploidy = sum(c(chrom_lengths$len[1:20]/1000*2) / sum(chrom_lengths$len / 1000))
  }
  
  # helper functions
  mytransform <- function(x,rho,phi) (phi*2^x-2*(1-rho))/rho
  geterrors <- function(rho,phi,meansSeg,weights, sds) {
    signal <- mytransform(meansSeg,rho,phi)
    mean(((round(signal)-signal)/sds)^2*weights/1000,na.rm=T)
  }
  
  get_purity = function(meansSeg, weights, normal_ploidy) {
    purs <- seq(0.05,1,.01)
    #ploidies <- seq(1.7,6,.02)
    ploidies = normal_ploidy
    errs <- matrix(NA, length(purs),length(ploidies))
    rownames(errs) <- purs
    colnames(errs) <- ploidies
    for(pp in 1:length(purs))
    {
      for(pl in 1:length(ploidies))
      {
        errs[pp,pl] <- geterrors(rho=purs[pp],phi=ploidies[pl],meansSeg,weights,(pl/pp)*0+1)
      }
    }
    mins <- arrayInd(which.min(errs),dim(errs))
    return(purs[mins[1]])
  }
  
  segmentation$purity = NA
  for (i in 1:nrow(segmentation)) {
    #TODO  calc purity for standard deviation to get uncertainty?
    segmentation$purity[i] <- get_purity(segmentation$value[i], 1, normal_ploidy)
  }
  
  # Experimental fitting procedure is disabled for now
  #segmentation$overall_purity = get_purity(segmentation$value, segmentation$len, normal_ploidy)
  segmentation$overall_purity = 1
  segmentation$assumed_ploidy = normal_ploidy
  segmentation$len = segmentation$len/1000000
  return(segmentation)
}
segmentation = get_purity_estimates(segmentation, sex)
write.table(segmentation, file=paste0("output/", samplename, "/", samplename,"_segmentation.txt"), quote=F, sep="\t", row.names=F)

max_value = max(abs(copynumber$logr_values))
if (max_value < 2) {
  max_y = 2
  signficance_bar_height = 0.05
} else if (max_value < 3) {
  max_y = 3
  signficance_bar_height = 0.05*1.5
} else {
  max_y = 5
  signficance_bar_height = 0.05*2.5
}


make_custom_plot = function(copynumber, segmentation, max_y, plot_title=NULL, plot_subtitle=NULL, signficance_bar_height=0.01) {
  background = data.frame(y=seq(max_y*(-1),max_y,1))
  p = ggplot(copynumber) +
    geom_hline(data=background, mapping=aes(yintercept=y), colour="black", alpha=0.3) +
    geom_point(mapping=aes(x=start, y=logr_values), alpha=0.5, size=0.9, colour="black") +
    geom_point(mapping=aes(x=start, y=segment_values), alpha=1, size=0.2, colour="red") +
    facet_grid(~chrom, scales="free_x", space = "free_x") +
    scale_x_continuous(expand=c(0, 0)) +
    ylim(-1*max_y,max_y) + ylab("log2 ratio") +
    theme_bw() + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       axis.text.y = element_text(colour="black",size=18,face="plain"),
                       axis.title.y = element_text(colour="black",size=20,face="plain"),
                       strip.text.x = element_text(colour="black",size=16,face="plain"),
                       plot.title = element_text(colour="black",size=36,face="plain",hjust = 0.5))
  if (any(segmentation$class=="loss")) {
    p = p + geom_rect(data=segmentation[segmentation$class=="loss",], mapping=aes(xmin=start, xmax=end, ymin=-(max_y-signficance_bar_height), ymax=-(max_y)), fill="blue")
  }
  if (any(segmentation$class=="gain")) {
    p = p + geom_rect(data=segmentation[segmentation$class=="gain",], mapping=aes(xmin=start, xmax=end, ymin=(max_y-signficance_bar_height), ymax=(max_y)), fill="red")
  }
  
  if (!is.null(plot_title) & !is.null(plot_subtitle)) {
    p = p + ggtitle(bquote(atop(.(plot_title), atop(.(plot_subtitle), ""))))
  } else if (!is.null(plot_title)) {
    p = p + ggtitle(plot_title)
  }
  return(p)
}


plot_title = dimnames(Biobase::assayDataElement(readCounts_tumour, "counts"))[[2]]
#plot_subtitle = paste0("Estimated purity: ", round(segmentation$overall_purity[1], 2), " Assumed ploidy: ", round(segmentation$assumed_ploidy, 2)) 
plot_subtitle = NULL
p = make_custom_plot(copynumber, segmentation, max_y=max_y, plot_title=plot_title, plot_subtitle=plot_subtitle, signficance_bar_height=signficance_bar_height)
png(file.path("output/", samplename, paste0(samplename, "_final_profile.png")), height=500, width=2000)
print(p)
dev.off()

################################################################################
#


segmentation$chrom = factor(segmentation$chrom, levels=gtools::mixedsort(unique(segmentation$chrom)))
segmentation$frac = segmentation$len / sum(segmentation$len)

breaks = seq(0, 1, 0.02)
segmentation$purity_group = cut(segmentation$purity, breaks)

purity_axis_breaks = levels(segmentation$purity_group)[(breaks*100) %% 20==0]
purity_axis_labels = breaks[(breaks*100) %% 20==0]


p_side_top = ggplot(segmentation) + aes(x=purity_group, y=frac) + geom_bar(position="stack", stat="identity") +
  scale_x_discrete(expand=c(0, 0.1), drop=FALSE, breaks=purity_axis_breaks, labels=purity_axis_labels) +
  theme_bw() + theme(axis.title.x=element_blank(),
                     axis.text.x=element_text(colour="black",size=14,face="plain", angle = 90, hjust = 1),
                     axis.ticks.x=element_blank(),
                     axis.text.y = element_text(colour="black",size=14,face="plain"),
                     axis.title.y = element_blank(),
                     strip.text.x = element_text(colour="black",size=16,face="plain"),
                     plot.title = element_text(colour="black",size=36,face="plain",hjust = 0.5),
                     plot.margin = margin(t=0.1, r=0.1, b=0.1, l=0.07, unit="cm"))
#unit(c(0.2,0.2,0.2,0.2), "cm"))


offset = 0.01
p_side_bottom = ggplot(segmentation) +
  #geom_hline(data=background, mapping=aes(yintercept=y), colour="black", alpha=0.3) +
  geom_rect(mapping=aes(xmin=start, xmax=end, ymin=purity-offset, ymax=purity+offset), alpha=0.5, size=0.9, colour="black", fill="black") +
  facet_grid(rows=vars(chrom), scales="free_y", space = "free_y", switch = "y") +
  scale_y_continuous(expand=c(0, 0), limits=c(0,1+offset), breaks=purity_axis_labels, labels=purity_axis_labels) +
  ylab("Fraction of cells") +
  theme_bw() + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     axis.text.x = element_text(colour="black",size=18,face="plain"),
                     axis.title.x = element_text(colour="black",size=20,face="plain"),
                     strip.text = element_text(colour="black",size=14,face="plain"),
                     plot.title = element_text(colour="black",size=36,face="plain",hjust = 0.5),
                     plot.margin = margin(t=0.1, r=0.1, b=0.1, l=0.1, unit="cm")) +
  coord_flip()

library(png)
library(grid)
rawdata = readPNG(file.path("output", samplename, paste0(samplename, "_rawdata.png")), native = FALSE)
png(file.path("output", samplename, paste0(samplename, "_combined.png")), height=1600, width=2000)
grid.arrange(arrangeGrob(rasterGrob(rawdata, interpolate = FALSE),
                         arrangeGrob(p_side_top, p_side_bottom, ncol=1, heights=c(1/5, 4/5)), ncol=2, widths=c(4/5, 1/5)),
             p, ncol=1, nrow=2, heights=c(0.7, 0.30))
dev.off()

save.image(file=paste0("output/", samplename, "/", samplename,"_output.RData"))
