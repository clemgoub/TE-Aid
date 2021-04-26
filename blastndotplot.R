
# self sequence db is being made by the shell script // or you need to make it to use in R
blastdotplot=function(query = NULL, db = NULL, getorf = NULL, blast = NULL, os = NULL){
  
  # run the selfblast
  bl=read.table(text=system(paste("blastn -query", query, "-db", db, "-evalue 0.05 -outfmt 6 -word_size 11 -gapopen 5 -gapextend 2 -reward 2 -penalty -3 | cut -f 1,7-10 | sed 's/#/-/g'"),
                            intern = TRUE)
                )
  # order from left to right
  bl=bl[order(bl$V2, decreasing = F),]

  # test if there are orf detected; store in orf if TRUE
  test<-try(read.table(as.character(getorf)), T)
   # test if there are TE prot detected; store in prot if TRUE
  test2<-try(read.table(as.character(blast)), T)
  #print(test)
  #print(class(test))
  if(class(test) == "data.frame"){
     orfs=read.table(as.character(getorf))
     }else{
     suppressWarnings(orfs$V1 <- as.data.frame(c(0,0))[,1])
   }

  #
  ## dot-plot (left)
  plot(x = 1, type = "n", xlim = c(0,bl$V3[1]),ylim = c(0,bl$V3[1]), col = "white",
       main = "TE consensus self dotplot (blastn)",
       ylab = paste(as.character(bl$V1[1]), "(bp)", sep = " "),
       xlab = paste(as.character(bl$V1[1]), "(bp)", sep = " ")
       )
    for(i in 1:length(bl$V1)){
      if(bl$V5[i] > bl$V4[i]){
        segments(x0 = bl$V2[i], x1 = bl$V3[i], y0 = bl$V4[i], y1 = bl$V5[i], col = "black", lwd = 1.5)
      } else {
        segments(x0 = bl$V2[i], x1 = bl$V3[i], y0 = bl$V4[i], y1 = bl$V5[i], col = "#009E73", lwd = 1.5)  
      }
        # if orientation
    } # for each segment end
  
  #
  ## Annotation graph (right)

  ## Arrows graph
  plot(x = 1, type = "n", xlim = c(0,bl$V3[1]),ylim = c(-max(length(orfs$V1),10),length(bl$V1)), col = "white",
       main = "TE consensus structure",
       xlab = paste(as.character(bl$V1[1]), "(bp)", sep = " "),
  )
  for(i in seq(1:length(bl$V1))){
    if(bl$V4[i] > bl$V2[i]){
        arrows(x0 = bl$V2[i], x1 = bl$V3[i], y0 = i, y1 = i, col = rainbow(length(bl$V1))[i], lwd = 3, length = 0.1)
        arrows(x0 = bl$V4[i], x1 = bl$V5[i], y0 = i, y1 = i, col = rainbow(length(bl$V1))[i], lwd = 3, length = 0.1)
    } # if to draw (filter)
  } # for each segment end

  ## orfs graphs
  if(class(test) != "data.frame"){
     text(paste("no orf >",os," bp detected", sep=""), x=bl$V3[1]/2, y=-5, cex = 2)
    } else {
      orfs=read.table(as.character(getorf))
      for(i in seq(1:length(orfs$V1))){
          rect(xleft = orfs$V2[i], xright = orfs$V3[i],
               ybottom = -i-0.05, ytop = -i+0.05, lwd = 1)
    }  # for each segment end
  } # if orf not null
  
  ## blastp graph
  if(class(test) != "data.frame"){
      prot=read.table(as.character(blast))
               
    } # If blast not empty
          ###### open blastp table
          ###### assign color for each class
          ###### draw colored rectangle same way as orf
          ###### print target name with text y = -i x = TE/2

    
} # function end


