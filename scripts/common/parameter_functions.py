def genrich_parameters(wildcards):
    if wildcards["suffix"] == "mm":
        return "-s 20"
