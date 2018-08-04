
ps 
ps.bact<-subset_taxa(ps, Kingdom=="Bacteria")
#trimming zeros if sum < 1000

ps.trim100.bact <- prune_taxa(taxa_sums(ps.bact) > 1000, ps.bact)
#agglomerate at family tax rank

ps.agg<-tax_glom(ps.trim100.bact, taxrank="Family", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
#prune baseline diet, sampling site (organ), and day

ps.bact.agg.2gr <- prune_samples(sample_data(ps.agg)$diet != "baseline", ps.agg) 
ps.bact.agg.2gr.col <- prune_samples(sample_data(ps.bact.agg.2gr)$organ == "CO", ps.bact.agg.2gr)
ps.bact.agg.2gr.col <- prune_samples(sample_data(ps.bact.agg.2gr.col)$day == "7", ps.bact.agg.2gr.col)
#plot the plot without ordination of taxa

plot_heatmap(ps.bact.agg.2gr.col, taxa.label = "Family", sample.order = "diet" )
#now we want to inform the abundance overview by arranging taxa on the heat-map according to their phylogenetics
#create a df that has all our Family taxa with their amplicon sequences

sequences <- colnames(ps.bact.agg.2gr.col@otu_table) 
at <- NULL
as.data.frame(at) 
taxanomic <- cbind(at,ps.bact.agg.2gr.col@tax_table) 
taxanomic <- as.data.frame(taxanomic) 
taxanomic$Genus <- rownames(taxanomic)

library(DECIPHER)
library(abind)
#constructor 
d <- taxanomic$Genus
d
z <- DNAStringSet(x=d, start=NA, end=NA, width=NA, use.names=TRUE)
z
names(z) <- rownames(taxanomic) # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(z), anchor=NA)

library(phangorn)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

#negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
#merge into the source Phyloseq object th etree was obtained from
library(phyloseq)
ps.mod.z <-merge_phyloseq(ps.bact.agg.2gr.col,fitGTR$tree)

ps.mod.z

taxa_names(ps) <- paste("SV",seq(nrow(tt.plus)))`
