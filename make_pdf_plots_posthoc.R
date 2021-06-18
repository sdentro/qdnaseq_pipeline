#
# This script creates a replot of the main result figure (and a few downsampled variants for easier manual editing) and it creates a plot of the CCF/purity values of called alterations
# The script can be supplied with a filtered set of calls after manual inspection and filtering. Each segment must then be annotated with the samplename and the inputfile can have all segments from all samples in the dataset concatenated
#

args = commandArgs(T)
filtered_calls_file = args[1]

if (!is.null(filtered_calls_file)) {
	final_calls = read.table(filtered_calls_file, header=T, stringsAsFactors=F)
}

library(ggplot2)
library(gridExtra)
setwd("output")
downsamplevalues = c(NA, 2, 4, 10)

for (samplename in list.files()) {
  print(samplename)
  
  setwd(samplename)
  load(paste0(samplename, "_output.RData"))
  copynumber_default = copynumber

  if (!is.null(filtered_calls_file)) {
    chrom_levels = levels(segmentation$chrom)
    segmentation = final_calls[final_calls$samplename==samplename,]
    segmentation$chrom = factor(segmentation$chromosome, levels=chrom_levels)
  }

  for (downsample in downsamplevalues) {
    
    outfilename = paste0(samplename, "_final_profile.pdf")
    
    if (!is.na(downsample)) {
      copynumber = copynumber_default[seq(1, nrow(copynumber_default), downsample),]
      outfilename = paste0(samplename, "_final_profile_sampleEvery", downsample, ".pdf")
    } else {
      copynumber = copynumber_default
    }

    plot_title = samplename    
    p = make_custom_plot(copynumber, segmentation, max_y=max_y, plot_title=plot_title, plot_subtitle=NULL, signficance_bar_height=signficance_bar_height)
    ggsave(p, file=file.path(outfilename), height=500/37.795, width=2000/37.795, dpi=300, useDingbats=FALSE, units="cm")
    
  }

  #=========================================================================================================================

    segmentation$chrom = factor(segmentation$chrom, levels=gtools::mixedsort(unique(segmentation$chrom)))
    segmentation$frac = segmentation$len / sum(segmentation$len)
    segmentation$purity[segmentation$classification=="normal"] = 0

    breaks = seq(0, 1, 0.02)
    segmentation$purity_group = cut(segmentation$purity, breaks)
    levels(segmentation$purity_group) = cut(breaks, breaks=breaks)
	
    # TODO: there is still a slight discrepancy in the breaks established here and the purity values in the other plot
    selected_breaks = which((breaks*100) %% 20==0)
    selected_breaks[length(selected_breaks)] = selected_breaks[length(selected_breaks)] - 2
    purity_axis_breaks = levels(segmentation$purity_group)[selected_breaks]
    #purity_axis_breaks = levels(segmentation$purity_group)[(breaks*100) %% 20==0]
    purity_axis_labels = breaks[(breaks*100) %% 20==0]


#print("all breaks")
#print(breaks)
#print("breaks")
#    print(purity_axis_breaks)
#print("labels")
#    print(purity_axis_labels) 

    mycolours = c("red", "blue", "white")
    names(mycolours) = c("gain", "loss", "normal")

    p_side_top = ggplot(segmentation) + aes(x=purity_group, y=frac, fill=classification) + geom_bar(position="stack", stat="identity", colour="black") +
      scale_x_discrete(expand=c(0, 0.1), drop=FALSE, breaks=purity_axis_breaks, labels=purity_axis_labels) +
      scale_fill_manual(values=mycolours) +
      ylab("Fraction of alterations") +
      theme_bw() + theme(axis.title.x=element_text(colour="black",size=14,face="plain"),
                         axis.text.x=element_text(colour="black",size=12,face="plain"),
                         axis.ticks.x=element_blank(),
                         axis.text.y = element_text(colour="black",size=14,face="plain"),
                         axis.title.y = element_blank(),
                         strip.text.x = element_text(colour="black",size=16,face="plain"),
                         plot.title = element_text(colour="black",size=36,face="plain",hjust = 0.5),
                         plot.margin = margin(t=0.1, r=0.1, b=0.1, l=0.07, unit="cm"),
			 legend.position="top",
			 legend.margin=margin(l = 0.1, b = 0, unit='cm'),
			 legend.text = element_text(size=14)) +
      guides(fill=guide_legend(title="Classification", title.theme = element_text(face="bold", size=10))) +
    coord_flip()
    #unit(c(0.2,0.2,0.2,0.2), "cm"))


    offset = 0.015
    p_side_bottom = ggplot(segmentation) +
      #geom_hline(data=background, mapping=aes(yintercept=y), colour="black", alpha=0.3) +
      geom_point(mapping=aes(x=start, y=0.001), size=0.00000000001, colour="white") + geom_point(mapping=aes(x=end, y=0.001), size=0.00000000001, colour="white") +
      geom_rect(mapping=aes(xmin=start, xmax=end, ymin=purity-offset, ymax=purity+offset, fill=classification), size=0.4, colour="black") +
      facet_grid(cols=vars(chrom), scales="free_x", space = "free_x") +
      scale_y_continuous(expand=c(0, 0), limits=c(0-2*offset,1+2*offset), breaks=purity_axis_labels, labels=purity_axis_labels) +
      scale_fill_manual(values=mycolours) +
      ylab("Fraction of cells") + xlab("Chromosome / Position") +
       theme_bw() + theme(axis.title.y=element_text(colour="black",size=14,face="plain"),
                         axis.text.y=element_text(colour="black",size=12,face="plain"),
                         axis.ticks.x=element_blank(),
			 axis.text.x = element_blank(),
                         axis.title.x = element_text(colour="black",size=14,face="plain"),
                         strip.text = element_text(colour="black",size=14,face="plain"),
                         plot.title = element_text(colour="black",size=36,face="plain",hjust = 0.5),
                         plot.margin = margin(t=0.4, r=0.1, b=0.55, l=0.1, unit="cm"),
			 panel.spacing = unit(0.01, "cm"),
			 legend.position="none")

    outfilename = paste0(samplename, "_ccf.pdf")
    # Disabled for now due to axis discrepancy between the two plots, see TODO above
    #p = arrangeGrob(p_side_bottom, p_side_top, ncol=2, widths=c(4/5, 1/5))
    p = p_side_bottom
    ggsave(p, file=file.path(outfilename), height=(500/37.795)/2, width=(2000/37.795)/2, dpi=300, useDingbats=FALSE, units="cm")


  setwd("../")
}
