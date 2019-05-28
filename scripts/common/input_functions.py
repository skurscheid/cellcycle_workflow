def chip_replicates(wildcards):
    fn = []
    for i in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"][wildcards["cycle"]]["ChIP"][wildcards["sample"]]:
        fn.append("/".join([wildcards["assayType"],
                            wildcards["project"],
                            wildcards["runID"],
                            "picardTools/MarkDuplicates",
                            i + wildcards["suffix"] + ".bam"]))
    return(fn)

def input_replicates(wildcards):
    fn = []
    for i in config["samples"]["ChIP-Seq"]["conditions"]["N08851_SK_LR1807201_SEQ"][wildcards["cycle"]]["Input"].values():
        for j in i:
            fn.append("/".join([wildcards["assayType"],
                                wildcards["project"],
                                wildcards["runID"],
                                "picardTools/MarkDuplicates",
                                j + wildcards["suffix"] + ".bam"]))
    return(fn)
