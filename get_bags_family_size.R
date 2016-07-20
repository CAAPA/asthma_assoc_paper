ped <- read.csv("data/input/fixed_pedigrees.txt", head=T, stringsAsFactors = F)[,c(1,2)]
adpc.sample.ids <- substring(read.table("data/input/chr22_bdos.dose.tfam", stringsAsFactors = F)[,1], 1, 16)
id.map <- read.delim("data/input/manifest_master.txt", head=T, stringsAsFactors = F)[,c(4,16)]
id.map <- id.map[id.map$Institute.Sample.Label %in% adpc.sample.ids,]
id.map <- merge(id.map, ped, by.x="Individual.ID", by.y="PATIENT", all.x=T)

length(unique(id.map$FAMILY))
mean(table(id.map$FAMILY))
