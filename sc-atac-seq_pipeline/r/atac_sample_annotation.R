#
# atac_sample_annotation.R
#
#


switch(cancer_type,
        "male-bc"={
                sampleAnnot <- data.frame(Sample = c("446B7L", "4CC61L"),
                          Color = sampleColors[1:2],
                          Cancer = c("male-bc", "male-bc"),
                          Histology = c("Infiltrating duct carcinoma", "Infiltrating duct carcinoma"),
                          BMI = c(26.0, 26.6),
                          Age = c(74, 65),
                          Race = c("AA", "AA"),
                          Stage = c("T1c", "T1c"),
                          Site = c("Breast", "Breast"),
                          Type = c("ER+ breast", "ER+ breast"))

        },
        "matt-ov"={
                sampleAnnot <- data.frame(Sample = c("3BAE2L", "3E5CFL"),
                          Color = sampleColors[1:2],
                          Cancer = c("ovarian", "ovarian"),
                          Histology = c("serous", "serous"),
                          BMI = c(22.13, 22.37),
                          Age = c(61, 59),
                          Race = c("CAU", "AS"),
                          Stage = c("IIB", "IIIC"),
                          Site = c("Ovary", "Ovary"),
                          Type = c("Ovarian", "Ovarian"))

        },
        {
                stop(sprintf("no reference data defined for %s", cancer_type))
        }
) # switch




