setwd("C:/Users/adolf/Desktop/LIB/Tesis")
rs_weight <- read.table("20k_rs.txt", header=T)
rs_interactors <- read.table("interacciones_reducido.txt", header=T)


final <- merge(x = rs_weight, y = rs_interactors, by.x = "sid", by.y = "Rs")
final_bc <- merge(x = rs_weight, y = rs_interactors, by.x = "sid", by.y = "Rs")

    
final$ldpred_beta <- final$ldpred_beta*(1+(final$NÂ._Interactores/608))

final <- final[c(2,3,1,4,5,6,7)]
final_bc <- final_bc[c(2,3,1,4,5,6,7)]

data <- data[c(1,3,2)]


write.table(final_bc, sep = "\t",
            file = "C:/Users/adolf/Desktop/LIB/Tesis/EUR.weight_LDpred_p1_final_bc.txt", 
            row.names = F, quote = F)

write.table(final, sep = "\t",
            file = "C:/Users/adolf/Desktop/LIB/Tesis/EUR.weight_LDpred_p1_final.txt", 
            row.names = F, quote = F)
