

out <- js.ms

#----- Model fit

library(plyr)

allchains2 <- rbind(as.matrix(out$samples[[1]]), 
                    as.matrix(out$samples[[2]]), 
                    as.matrix(out$samples[[3]]))

z.act <- grep("zzz.fit", colnames(allchains2))[1]
z.new <- grep("zzz.fit", colnames(allchains2))[2]

p2 <- round(mean(allchains2[,z.act]>allchains2[,z.new]),2)
m3 <- round_any(min(allchains2[,z.new],allchains2[,z.act]), 10, f = floor)
m4 <- round_any(max(allchains2[,z.new],allchains2[,z.act]), 10, f = ceiling)
plot(allchains2[,z.act], allchains2[,z.new], xlab = expression(T^{obs}), ylab=expression(T^{rep}), cex.lab=1, cex.axis=1, xlim = c(m3,m4), ylim = c(m3,m4), las=1, main = "Modified Pearson Residuals")
abline(0, 1, lwd=2); mtext(p2, side = 3, line = -2, at=(m4-m3)*.1 + m3, cex = 2)








