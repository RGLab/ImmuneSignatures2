library(ImmuneSignatures2)
library(ImmuneSpaceR)
library(Rlabkey)
library(data.table)

studies <- c("SDY56, SDY61, SDY63, SDY67, SDY80, SDY180, SDY212, SDY224, SDY269,
             SDY270, SDY400, SDY404, SDY520, SDY640, SDY984, SDY1119, SDY1260, SDY1264,
             SDY1276, SDY1289, SDY1291, SDY1293, SDY1294, SDY1325, SDY1364, SDY1368, SDY1370,
             SDY1373")
studies <- strsplit(studies, ", ")[[1]]

con <- CreateConnection("", onTest = TRUE)

inputSmpls <- getTable(con, "assay.ExpressionMatrix.matrix", "inputSamples_computed")
colnames(inputSmpls) <- gsub("^biosample_","", colnames(inputSmpls))
inputSmpls$study <- gsub("^SUB\\d{6,8}\\.", "SDY", inputSmpls$participantid)
inputSmpls <- inputSmpls[ inputSmpls$study %in% studies, ]

demographics <- con$getDataset("demographics")
demographics <- addStudy(demographics)
demographics <- imputeAge(demographics)

inputSmpls$age_imputed <- demographics$age_imputed[ match(inputSmpls$participantid, demographics$participant_id)]
inputSmpls$age_bracket <- ceiling(inputSmpls$age_imputed/10) * 10

inputSmpls <- data.table(inputSmpls)
res <- inputSmpls[, list(count = .N), by = c("age_bracket", "study")]

# stacked bar plot
library(ggplot2)
ggplot(res, aes(fill=study, y=count, x=age_bracket)) +
  geom_bar(position="stack", stat="identity")
