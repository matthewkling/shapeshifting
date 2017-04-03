

# "Seasons Shifting Shape" visualization
# Matthew Kling, April 2017
# mattkling@gmail.com

# prior to running script, download DAYMET's "multiple extractor script"
# from https://daymet.ornl.gov/web_services.html into R working directory.


library(data.table)
library(dplyr)
library(tidyr)
library(jsonlite)
library(grid)
library(gridExtra)
library(stringr)
library(ggplot2)
library(ggmap)



# query daymet API for climate data
get_ingredients <- function(cities){

      # lat-lon lookup
      points <- geocode(cities, source="google")
      
      # create a lat-lon file for daymet extractor
      cities <- substr(cities, 1, regexpr(",", cities)-1)
      dir.create("daymet_data")
      files <- paste0("daymet_data/", cities, ".csv")
      write.table(cbind(files, points[,2], points[,1]), 
                  "city_coordinates.txt",
                  col.names=F, row.names=F, sep=", ", quote=F)
      
      # call daymet script to query API
      system("java -Xms512m -Xmx1024m -Dhttps.protocols=TLSv1.1,TLSv1.2 -jar daymet_multiple_extraction.jar city_coordinates.txt")
}

# generate individual charts for single locations
make_sausage <- function(file, window=31, cutoff=1995, plot=T){
      
      # load data
      d <- fread(file) %>%
            as.data.frame() %>%
            tbl_df()
      names(d) <- str_trim(substr(names(d), 1, 4))
      
      # distill data. this could be cleaned up.
      m <- as.data.frame(d) %>%
            filter(year<cutoff) %>%
            group_by(yday) %>%
            summarize(ppt=mean(prcp, na.rm=T),
                      tmin=mean(tmin, na.rm=T),
                      tmax=mean(tmax, na.rm=T)) %>%
            mutate(tmean=(tmin+tmax)/2,
                   ppt=as.vector(stats::filter(ppt, rep(1, window), 
                                               method="convolution",
                                               sides=2, circular=T)),
                   tmean=as.vector(stats::filter(tmean, rep(1, window), 
                                                 method="convolution", 
                                                 sides=2, circular=T)),
                   
                   ppt=as.vector(stats::filter(ppt, rep(1, window), 
                                               method="convolution",
                                               sides=2, circular=T)),
                   tmean=as.vector(stats::filter(tmean, rep(1, window), 
                                                 method="convolution", 
                                                 sides=2, circular=T)),
                   
                   period="full")
      
      m2 <- as.data.frame(d) %>%
            filter(year>=cutoff) %>%
            group_by(yday) %>%
            summarize(ppt=mean(prcp, na.rm=T),
                      tmin=mean(tmin, na.rm=T),
                      tmax=mean(tmax, na.rm=T)) %>%
            mutate(tmean=(tmin+tmax)/2,
                   
                   ppt=as.vector(stats::filter(ppt, rep(1, window), 
                                               method="convolution",
                                               sides=2, circular=T)),
                   tmean=as.vector(stats::filter(tmean, rep(1, window), 
                                                 method="convolution", 
                                                 sides=2, circular=T)),
                   
                   ppt=as.vector(stats::filter(ppt, rep(1, window), 
                                               method="convolution",
                                               sides=2, circular=T)),
                   tmean=as.vector(stats::filter(tmean, rep(1, window), 
                                                 method="convolution", 
                                                 sides=2, circular=T)),
                   period="21C")
      
      m <- rbind(m, m2)
      m <- tbl_df(m)
      
      s <- m[m$yday %in% (c(0:365)*1),] %>%
            dplyr::select(yday, ppt, tmean, period) %>%
            gather(var, value, ppt, tmean) %>%
            unite(stat, var, period) %>%
            spread(stat, value)
      
      m <- rbind(m, filter(m, yday==1)) # close the loop
      
      p <- ggplot() + 
            geom_path(data=filter(m, period=="full"), 
                      aes(tmean, ppt, 
                          order=yday, 
                          group=period,
                          color=yday),
                      size=.1, alpha=.5) +
            geom_segment(data=s, 
                         aes(x=tmean_full, y=ppt_full, 
                             xend=tmean_21C, yend=ppt_21C,
                             group=yday, color=yday),
                         size=.1) +
            geom_path(data=filter(m, period=="21C"), 
                      aes(tmean, ppt, 
                          order=yday, color=yday, 
                          group=period, alpha=period),
                      size=1, lineend="round") +
            annotate(geom="text", 
                     label=sub(".csv", "", basename(file)), 
                     size=7, hjust=.5, color="gray40",
                     x=mean(range(c(s$tmean_full, s$tmean_21C))),
                     y=max(c(s$ppt_full, s$ppt_21C)) + 
                           diff(range(c(s$ppt_full, s$ppt_21C)))/10) +
            scale_color_gradientn(colours=c("dodgerblue", "cyan", "green", 
                                            "yellow", "red", "magenta", "dodgerblue")) +
            scale_alpha_manual(values=c(1, .4)) +
            scale_x_continuous(expand=c(.05,0)) +
            scale_y_continuous(expand=c(.05,0)) +
            theme_minimal() +
            theme(panel.background=element_rect(fill="black", color="black"),
                  panel.grid=element_blank(),
                  axis.text=element_blank(),
                  axis.title=element_blank(),
                  plot.background=element_rect(fill="black", color="black"),
                  plot.title=element_text(color="gray40", size=20),
                  legend.position="none") +
            labs(x="temperature",
                 y="precipitation")
      
      if(plot){
            dir.create("subplots")
            ggsave(paste0("subplots/", sub("csv", "pdf", basename(file))),
                      p, width=6, height=6, units="in")
      }
      
      return(p)
      
      
}

# compile final visualization
package_sausages <- function(plots, outfile){
      
      ### additional chart components
      
      pad <- ggplot() + theme(panel.background=element_rect(fill="black", color="black"),
                              panel.grid=element_blank(),
                              axis.text=element_blank(), axis.ticks=element_blank(),
                              plot.background=element_rect(fill="black", color="black"))
      
      title <- pad +
            annotate(geom="text", 
                     label="THE SHIFTING\nSHAPE OF\nTHE SEASONS", size=20,
                     hjust=.5, color="gray40",
                     x=0, y=0, fontface="bold", lineheight=.75)
      
      subtitle <- pad +
            annotate(geom="text", 
                     label="annual cycles of daily\ntemperature & precipitation\nover space & time", 
                     size=12, hjust=.5, color="gray40",
                     x=0, y=0, lineheight=.75)
      
      credits <- pad +
            annotate(geom="text", 
                     label="Visualization by Matthew Moore Kling.\nCode at github.com/matthewkling/shapeshift.\nDaily climate data from DAYMET.", 
                     size=7, hjust=.5, color="gray20",
                     x=0, y=0, lineheight=1)
      
      
      # expository legend
      ld <- data.frame(day = 1:365) %>%
            mutate(angle = day/365*2*pi,
                   radius = 10,
                   y = sin(-angle + pi/2) * radius, # orient like clock and convert to xy
                   x = cos(-angle + pi/2) * radius)
      ld <- full_join(ld, mutate(ld, 
                                 x2=x+2, # offset time periods
                                 y2=y+5))
      ld <- rbind(ld, ld[1,]) # close the loop
      ml <- data.frame(letter = c("J", "F", "M", "A", "M", "J",
                                  "J", "A", "S", "O", "N", "D"),
                       month = (1:12) - .5, # center lables on month, unlike a clock
                       stringsAsFactors = F) %>%
            mutate(angle = month/12*2*pi,
                   radius = 9,
                   y = sin(-angle + pi/2) * radius,
                   x = cos(-angle + pi/2) * radius,
                   x = x + 2,
                   y = y + 5)
      
      expo <- ggplot() +
            annotate(geom="text", 
                     x=c(-6, -2)+9, 
                     y=c(-9.5,-3)+17, 
                     label=c("1980-\n         1999", 
                             "2000-\n         2015"),
                     fontface=c("plain", "bold"),
                     alpha=c(.6, 1),
                     size=c(6,8),
                     hjust=.5, vjust=0,
                     color="dodgerblue", lineheight=.75) +
            geom_path(data=ld, 
                      aes(x, y, color=angle),
                      size=.1, alpha=.5, lineend="round") +
            geom_segment(data=distinct(ld),
                         aes(x, y, xend=x2, yend=y2, color=angle),
                         size=.1, lineend="round") +
            geom_path(data=ld, 
                      aes(x2, y2, color=angle),
                      size=1, lineend="round") +
            geom_text(data=ml,
                      aes(x, y, color=angle, label=letter),
                      size=8) +
            scale_color_gradientn(colours=c("dodgerblue", "cyan", "green", 
                                            "yellow", "red", "magenta", "dodgerblue")) +
            annotate(geom="text", 
                     x=c(-0, -12), 
                     y=c(-12, 0), 
                     angle=c(0,90),
                     label=c("-    temperature    +",
                             "-    precipitation    +"),
                     fontface="bold",
                     hjust=.5,
                     color="gray40", size=8, lineheight=.75) +
            scale_alpha_manual(values=c(1, .4)) +
            scale_x_continuous(expand=c(.05,0)) + # use to control dims independent of grob layout
            scale_y_continuous(expand=c(.1,0)) +
            theme_minimal() +
            theme(panel.background=element_rect(fill="black", color="black"),
                  panel.grid=element_blank(),
                  axis.text=element_blank(),
                  axis.title=element_blank(),
                  plot.background=element_rect(fill="black", color="black"),
                  plot.title=element_text(color="gray40", size=20),
                  legend.position="none") +
            labs(x="temperature",
                 y="precipitation")
      
      
      ### page layout
      
      # main content
      lo <- read.csv("e:/shapeshift/city_names.csv", header=F) %>%
            as.matrix() %>% t() %>% as.vector()
      lo <- substr(lo, 1, regexpr(",", lo)-1)
      lo <- match(lo, sort(lo))
      mg <- arrangeGrob(grobs=plots[lo], nrow=6, clip=T)
      mg <- arrangeGrob(pad, mg, pad, 
                        heights=c(1,10,1), ncol=1)
      
      # metadata
      tg <- arrangeGrob(pad, title, subtitle, credits, pad, expo, pad, pad, 
                        heights=c(2, 1.3, 1, .9, 1.2, 3.4, 1.3,  1), ncol=1)
      
      # final arrangement
      g <- arrangeGrob(pad, tg, pad, mg, pad, 
                       widths=c(.75,2,.5,10,.75), nrow=1)
      
      # export
      pdf(outfile, 
          width=50, height=30,
          bg="black", fg="black")
      grid.draw(g)
      dev.off()
      
}



####### RUN ######

cities <- c("Seattle, WA", "Spokane, WA", "Kalispell, MT", "Bismarck, ND", "Minneapolis, MN", "Milwaukee, WI", "Detroit, MI", "Bangor, ME",
            "Portland, OR", "Bend, OR", "Bozeman, MT", "Rapid City, SD", "Omaha, NB", "Chicago, IL", "Pittsburgh, PA", "Boston, MA", 
            "Eureka, CA", "Boise, ID", "Idaho Falls, ID", "Casper, WH", "Kansas City, KS", "Indianapolis, IN", "Washington, DC", "New York, NY", 
            "San Francisco, CA", "Reno, NV", "Salt Lake City, UT", "Denver, CO", "Tulsa, OK", "Saint Louis, MO", "Charlotte, NC", "Virginia Beach, VA",
            "Los Angeles, CA", "Las Vegas, NV", "Flagstaff, AZ", "Albuquerque, NM", "Dallas, TX", "Birmingham, AL", "Atlanta, GA", "Charleston, SC",
            "San Diego, CA", "Phoenix, AZ", "Tucson, AZ", "El Paso, TX", "Houston, TX", "New Orleans, LA", "Tampa, FL", "Miami, FL")

#get_ingredients(cities) # run this once to get data

list.files("daymet_data", full.names=T) %>%
      lapply(make_sausage, window=31, cutoff=2000, plot=F) %>%
      package_sausages(outfile="shapeshifting.pdf")

