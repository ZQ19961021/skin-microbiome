 
library(vegan) 
library(ape)
library(ggplot2)
library(RColorBrewer)
library(ade4)
library(ggrepel)


inf="species0619.profile"
inp='mapping_file.txt'

species<-function(ii){
  ii<-as.character(ii)
  ii[intersect(grep('.*[|,.;]g_.*[|,.;]s_..*',ii),grep('[|,.;]s_$',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]g_.*[|,.;]s_..*',ii),grep('[|,.;]s_$',ii,invert=T))],function(x){gsub('.*[|,.;]s','s',x)}))
  ii[intersect(grep('.*[|,.;]f_.*[|,.;]g_..*',ii),grep('g_[|,.;]s_',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]f_.*[|,.;]g_..*',ii),grep('g_[|,.;]s_',ii,invert=T))],function(x){gsub('.*[|,.;]g','g',x)}))
  ii[intersect(grep('.*[|,.;]o_.*[|,.;]f_..*',ii),grep('f_[|,.;]g_',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]o_.*[|,.;]f_..*',ii),grep('f_[|,.;]g_',ii,invert=T))],function(x){gsub('.*[|,.;]f','f',x)}))
  ii[intersect(grep('.*[|,.;]c_.*[|,.;]o_..*',ii),grep('o_[|,.;]f_',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]c_.*[|,.;]o_..*',ii),grep('o_[|,.;]f_',ii,invert=T))],function(x){gsub('.*[|,.;]o','o',x)}))
  ii[intersect(grep('.*[|,.;]p_.*[|,.;]c_..*',ii),grep('c_[|,.;]o_',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('.*[|,.;]p_.*[|,.;]c_..*',ii),grep('c_[|,.;]o_',ii,invert=T))],function(x){gsub('.*[|,.;]c','c',x)}))
  ii[intersect(grep('k_.*[|,.;]p_..*',ii),grep('k_[|,.;]p_',ii,invert=T))]<-as.character(lapply(ii[intersect(grep('k_.*[|,.;]p_..*',ii),grep('k_[|,.;]p_',ii,invert=T))],function(x){gsub('.*[|,.;]p','p',x)}))
  return(ii)
}

prof <- read.csv(file=inf,sep="\t",row.names=1,check.names = F,stringsAsFactors = F)  
#prof = prof[grep("Viruse",row.names(prof)),]
#row.names(prof) = species(row.names(prof))
map <- read.csv(file=inp,sep="\t",check.names = F,stringsAsFactors = F)
map = map[grep("DP",map[,1]),]
#map = map[map$Type=="Fh",]
#map = map[map$Gender=="Female",]
prof = prof[,colnames(prof)%in%map[,1]]
row.names(map) = map[,1]
map = map[colnames(prof),]
prof <- t(prof)
prof <- prof[,colSums(prof)!=0]
adonis(prof~map$age_type)
prof <- sqrt(prof)

#dbRDA========================================================================================
ord <- capscale(prof~map$Age,map,distance = "bray") 
sp = scores(ord,choices = 1:2,display = "sites")
s.class(sp,factor(map[,5]),col = brewer.pal(4,"Set1"))

pFIG <- ordisurf(ord,map[,'Age'], bubble = 1) 

extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz)
}
contour.vals <- extract.xyz(obj = pFIG)

p<-plot(ord)
site.scores <- as.data.frame(scores(p, "site")) 
species.scores <- as.data.frame(scores(p, "species"))
site.scores$lab <- row.names(site.scores)
site.scores <- cbind(site.scores,map)[-4]
site.scores$z <- NA 

x1 <- min(site.scores[,1]) - 0.3
y1 <- min(site.scores[,2]) - 0.3
x2 <- max(site.scores[,1]) + 0.3
y2 <- max(site.scores[,2]) + 0.3
inlength =1
bb <- head(species.scores[order(sqrt((species.scores[,1])^2+(species.scores[,2])^2),decreasing=T),],n=20L)[1:20,]
cutoff <- (x2-0.3) / abs(bb[1,1]) * inlength
d2 <- data.frame(X=bb[1:dim(bb)[1],1]*cutoff, Y=bb[1:dim(bb)[1],2]*cutoff, LAB=rownames(bb)[1:dim(bb)[1]])
d2$z = NA
color = c(brewer.pal(4,"Set1"))
site.scores$color = color[1]
site.scores$color[site.scores$age_type=="old"] = color[2]
site.scores$color[site.scores$age_type=="older"] = color[3]
#site.scores$color[site.scores$group=="Longevity"] = color[4]
ggplot(data = contour.vals, aes(x, y, z=z)) +
  stat_contour(aes(colour = ..level..))+
  scale_color_continuous(low='grey80',high = '#272EAF')+
  geom_point(data = site.scores,size=2,
             aes(x=CAP1,y=MDS1,z=z,shape=site.scores$Type),colour = site.scores$color)+
  #geom_text_repel(data=d2, aes(x=X, y=Y, label=LAB),
  #               family="Helvetica", fontface="italic", size=3, check_overlap=TRUE) +
  #geom_segment(data=d2, aes(x=0, y=0, xend=X, yend=Y),
  #             arrow = arrow(length = unit(0.3, "cm")), size=0.8, alpha=0.5)+
  theme_bw()+
  coord_cartesian(ylim = c(-2, 2))+
  labs(title="", x="CAP1", y="MDS1")+
  theme(
    panel.grid = element_blank()
  )



