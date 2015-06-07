dir.create("animation")
setwd("animation")

# example 1: simple animated countdown from 10 to "GO!".
png(file="example%02d.png", width=893, height=721)
vec = 1:1000
for(i in 1:1000){
  values = test_cont.liste[[1]][i,]
  header = sprintf("SD = %3.2f, Q0 = %3.2f, Kappa = %3.2f, Vs30 = %3.2f, MAG = %3.2f, Dist = %3.2f, PGA = %3.2f\n", 
                   values$SD, values$Q0, values$kappa, values$Vs30, values$MAG, values$DIST, values$PGA)
  
  mp = barplot(m[2,i,], 
               main = header,
               ylab = "Probability",
               xlab = "log(PGA)",
               axes = F, 
               ylim = c(0,1.05), plot = F)
  
  axis(1, at = mp, labels = midpoints)
  axis(2, at = seq(0,1,0.1), labels = seq(0,1,0.1), las = 2)
  box()
  
  vec[i] = max(midpoints[which(max(m[2,i,])==m[2,i,])]) 
  
  #Sys.sleep(0.5)
  #readline(prompt="")
}
dev.off()

# convert the .png files to one .gif file using ImageMagick. 
# The system() function executes the command as if it was done
# in the terminal. the -delay flag sets the time between showing
# the frames, i.e. the speed of the animation.
shell(paste("cd",getwd(), '&& "C:\\Program Files\\ImageMagick-6.9.1-Q16\\convert.exe" -delay 80 *.png example_1.gif',sep=" "), translate = T)

# to not leave the directory with the single jpeg files
# I remove them.
file.remove(list.files(pattern=".png"))
