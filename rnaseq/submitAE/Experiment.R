all.out <- list()
all.out[["Comment[ArrayExpressAccession]"]] <- "E-MTAB-5308"
all.out[["MAGE-TAB Version"]] <- "1.1"
all.out[["Investigation Title"]] <- "Quantifying the transcriptional effects of lncRNA depletion in HeLa cells with RNA-seq"
all.out[["Comment[Submitted Name]"]] <- "Quantifying the transcriptional effects of lncRNA depletion in HeLa cells with RNA-seq"
all.out[["Experiment Description"]] <- "Long noncoding RNAs (lncRNAs) are a major transcriptional output of the mammalian genome, yet their functions are poorly characterized. In this experiment, we used RNA interference (RNAi) and locked nucleic acid antisense oligonucleotides (LNAs) to deplete the lncRNAs LINC00899 and C1QTNF1-AS1 in HeLa cells. Both of these lncRNAs are of interest as their depletion results in mitotic delay, suggesting a role in mitotic progression. After depletion of each lncRNA with RNAi or LNAs, we performed high-throughput RNA sequencing and tested for differential expression compared to negative control samples (using control siRNAs or negative control LNAs) to identify downstream regulatory targets of each lncRNA. We generated 3-4 replicates for each condition, spread across multiple batches and experiments. Samples with the same batch number were treated at roughly the same time (within the same month), while the experiment number is nested within batch and refers to cells treated on the same day. Samples were pooled and sequenced across multiple lanes, yielding 8-9 technical replicates per sample."

all.out[["Experimental Design"]] <- c("genotype design",
                                      "compound treatment design",
                                      "cellular modification design")
all.out[["Experimental Design Term Source REF"]] <- c("EFO", 
                                                      "EFO",
                                                      "EFO")
all.out[["Experimental Design Term Accession Number"]] <- c("EFO_0001748",
                                                            "EFO_0001755",
                                                            "EFO_0004666")

all.out[["Experimental Factor Name"]] <- c("compound", 
                                           "genotype",
                                           "loss of function method")
all.out[["Experimental Factor Type"]] <- all.out[["Experimental Factor Name"]]
all.out[["Experimental Factor Term Source REF"]] <- ""
all.out[["Experimental Factor Term Accession Number"]] <- ""

all.out[["Person Last Name"]] <- "Lun"
all.out[["Person First Name"]] <- "Aaron"
all.out[["Person Mid Initials"]] <- "TL"
all.out[["Person Email"]] <- "aaron.lun@cruk.cam.ac.uk"
all.out[["Person Phone"]] <- ""
all.out[["Person Fax"]] <- ""      
all.out[["Person Address"]] <- "University of Cambridge Li Ka Shing Centre Robinson Way Cambridge CB2 0RE United Kingdom"
all.out[["Person Affiliation"]] <- "Cancer Research UK Cambridge Institute"
all.out[["Person Roles"]] <- "submitter"

all.out[["Protocol Name"]] <- c("P-MTAB-53234", 
                                "P-MTAB-53235",
                                "P-MTAB-53236",
                                "P-MTAB-53237",
                                "P-MTAB-53238",
                                "P-MTAB-53239",
                                "P-MTAB-53240",
                                "P-MTAB-53241")
all.out[["Protocol Type"]] <- c("sample collection protocol",
                                "growth protocol",
                                "treatment protocol",
                                "nucleic acid extraction protocol",
                                "nucleic acid library construction protocol",
                                "nucleic acid sequencing protocol",
                                "high throughput sequence alignment protocol",
                                "conversion protocol")
all.out[["Protocol Term Source REF"]] <- c("EFO", 
                                           "EFO",   
                                           "EFO",
                                           "EFO",
                                           "EFO",
                                           "EFO",
                                           "EFO",
                                           "EFO")
all.out[["Protocol Term Accession Number"]] <- c("EFO_0005518",
                                                 "EFO_0003789",
                                                 "EFO_0003969",
                                                 "EFO_0002944",
                                                 "EFO_0004184",
                                                 "EFO_0004170",
                                                 "EFO_0004917",
                                                 "EFO_0005520")

all.out[["Protocol Description"]] <- c("HeLa cells were obtained from American Type Culture Collection.", 
                                       "HeLa cells were maintained in Dulbecco's modified Eagle's medium (Sigma Aldrich, D6429) supplemented with 10% fetal bovine serum (Thermo Fisher Scientific) and cultured at 37 degrees Celsius with with 5% CO2.",
                                       "HeLa cells were transfected with Lipofectamine RNAiMax reagent (Thermo Fischer Scientific) following the manufacturer's instructions, using siRNAs (Thermo Fischer Scientifc) and LNA Gapmers (Exiqon) at a final concentration of 50 nM and 25nm, respectively. All experiments were done 48 hr after transfection.",
                                       "RNA (1 ug) was extracted with the RNeasy Kit (QIAGEN) and treated with DNase I following the manufacturer's instructions. RNA quality was assessed using a Total RNA Nano chip with a 2100 Bioanalyzer instrument (Agilent).",
                                       "RNA-seq libraries were prepared from HeLA cells using TruSeq Stranded Total RNA Kit with Ribo-Zero Gold (Illumina, RS-122-2303). Library quality was assessed using a DNA1000 chip with a 2100 Bioanalyzer instrument (Agilent).",
                                       "Indexed libraries were PCR amplified and sequenced on multiple lanes of an Illumina Hiseq 2500 instrument to obtain 125 bp paired-end reads.",
                                       "Reads were aligned to the hg38 build of the human genome using subread v1.6.1 in paired-end RNA-seq mode. The number of read pairs mapped to the exonic regions of each gene was then counted for each library, using the featureCounts function in Rsubread v1.30.3 with Ensembl GRCh38 version 91. Only alignments with mapping quality scores above 10 and forming reversely stranded fragments were considered during counting.",
                                       "The QuantiTect Reverse Transcription Kit (QIAGEN) was used for cDNA synthesis including an additional step to eliminate genomic DNA contamination.")
all.out[["Protocol Hardware"]] <- c("BD FACSAria III cell sorter",
                                    "", 
                                    "",
                                    "2100 Bioanalyzer",
                                    "2100 Bioanalyzer",
                                    "Illumina Hiseq 2500, Illumina HiSeq 4000",
                                    "",
                                    "")
all.out[["Protocol Software"]] <- c("", 
                                    "",
                                    "",
                                    "",
                                    "",
                                    "",
                                    "(R)subread",
                                    "")
                                 
all.out[["Term Source Name"]] <- "EFO"
all.out[["Term Source File"]] <- "http://www.ebi.ac.uk/efo/"
all.out[["Term Source Version"]] <- ""
all.out[["Public Release Date"]] <- "2017-05-31"
all.out[["Comment[AEExperimentType]"]] <- "RNA-seq of coding RNA"
all.out[["Comment[SequenceDataURI]"]] <- "http://www.ebi.ac.uk/ena/data/view/ERR1751045-ERR1751697"
all.out[["Comment[SecondaryAccession]"]] <- "ERP020478"
all.out[["SDRF File"]] <- "E-MTAB-5308.sdrf.txt"

unlink("idf.tsv")
for (x in names(all.out)) {
    write(file="idf.tsv", paste0(c(x, all.out[[x]]), collapse="\t"), append=TRUE)
}


