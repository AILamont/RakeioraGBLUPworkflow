#workflow
rule rscript:
        input:
            Rfile = "GBLUPsnakemakeMay23.R",
            PED="../InputFiles/simped103id",
            VCF="../InputFiles/chr22_filteredCHB2thin.recode.vcf",
            PHENO="../InputFiles/simdatafull"
        output:
            csv = "../Results/snakedevout.csv"
        threads: 1
        shell:
            "Rscript {input.Rfile} GBLUP {input.PED} {input.VCF} {input.PHENO} NA > {output.csv}"