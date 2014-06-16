#########################
#
# Copyright (c) 2013, SK Woolf <bcnskaa@gmail.com>.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#  1. Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
# 
#  2. Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
#########################
source("process_DF.R");

# https://personality-project.org/r/r.graphics.html

## the following code and figure is adapted from the help file for pairs 
plot_scatter_mtx2 <- function(plot_data, plot_title=character(0), pdf_outfile=character(0), png_resolution=120, png_outfile=character(0), scatter_size=100) 
{
	library(graphics);
	data_n <- dim(plot_data)[2];
	
	if(length(pdf_outfile) != 0)
	{
		print(paste("Plotting ", plot_title, " to file ", pdf_outfile, sep=""));
		pdf_width <- ((data_n * scatter_size) + 100) / 75;
		pdf_height <- ((data_n * scatter_size) + 100) / 75;	
		pdf(pdf_outfile, width=pdf_width, height=pdf_height);
	} else if(length(png_outfile) != 0) {
		png_width <- ((data_n * scatter_size) + 100) / 75;
		png_height <- ((data_n * scatter_size) + 100) / 75;	
		
		#print(paste("PNG image size=", png_width, "x", png_height, sep=""));
		png(png_outfile, width=png_width, height=png_height, pointsize=4, res=png_resolution, unit="in");
		#png(png_outfile);
	}

	
	print("Plotting data...");
	par(pch=20);
	par(cex=0.5);
	pairs(plot_data, lower.panel=panel.xylm, upper.panel=panel.cor, main=plot_title);
	#dev.off();
	
	if(length(pdf_outfile) != 0 | length(png_outfile) != 0)
	{
		dev.off();
	}
}

plot_scatter_mtx <- function(bin_list, plot_title=character(0), pdf_outfile=character(0), png_resolution=120, png_outfile=character(0), scatter_size=100, trim_plot_data=T) 
{	
	
	data_n <- length(bin_list);
	labels <- names(bin_list);
	row_n <- as.numeric(unique(sapply(bin_list, function(v) length(v))))
	
	# plot parameters
	if(length(plot_title) == 0)
		plot_title <- paste("Correlation of ", data_n, " datasets", sep="");
	
	print(paste("data_n=", data_n, " and row_n=", row_n, sep=""))
	if(length(row_n) > 1)
	{
		print("Bin number is not consistent, unable to plot the data.");
		return;
	}
	

	#rownames(plot_data) <- rownames(bin_list);
	plot_data <- create_plot_data(bin_list);
	
	if(trim_plot_data == T)
	{
		print("Trimming plot data...");
		plot_data <- trim_plot_data(plot_data);	
	}
	
	plot_scatter_mtx2(plot_data, plot_title=plot_title, pdf_outfile=pdf_outfile, png_resolution=png_resolution, png_outfile=png_outfile, scatter_size=scatter_size);
}




## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
## first create a function (panel.cor)
## then create a scatterplot
## Usage:
## pairs(plot_data2, lower.panel=panel.xylm, upper.panel=panel.cor, main="Main Title");
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
	rgb_palette_max <- 1;
	rgb_palette_min <- -1;
	rgb_palette_step <- 2000;
	rgb_palette_leading <- 1 / (rgb_palette_step * rgb_palette_step);
	rgb.palette <- colorRampPalette(c("red", "white", "blue"), space = "rgb")(rgb_palette_step);
	
	usr <- par("usr"); on.exit(par(usr));
	par(usr = c(0, 1, 0, 1));
	r = (cor(x, y));
	txt <- format(c(r, 0.123456789), digits=digits)[1];
	txt <- paste(prefix, txt, sep="");
	if(missing(cex.cor)) cex <- 1.0 / strwidth(txt);
	#text(0.5, 0.5, txt, cex = cex * abs(r))
	
	#print(paste("sum(x)=", sum(x), ", sum(y)=", sum(y), sep=""));
	#print(paste("cor(", colnames(x), ", ", colnames(y), ")=", r, sep=""));
	
	if(r == rgb_palette_max) {
		r_col <- rgb.palette[rgb_palette_step];
	} else {
		r_col <- rgb.palette[ceiling(((r - rgb_palette_min) + rgb_palette_leading) / (1 / (rgb_palette_step / 2)))];
	}
	
	#print(paste("Text=", txt, sep=""))
	rect(0,0,1,1,col=r_col);
	text(0.5, 0.5, txt, cex = 2)
}



panel.xylm<-function (x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 0.5, col.lm = "green", lwd=par("lwd"), ...)
{
	#print(paste("Plotting xylm..."));
	points(x, y, pch = 20, col = col, bg = bg, cex = cex);
	
	ok <- is.finite(x) & is.finite(y)
	if (any(ok)) {
		abline(lm(y~x,subset=ok), col = col.lm, ...);
		#abline(lowess(y~x,subset=ok), col = col.lm, ...);
	}
}

panel.xylm_smooth <-function (x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 0.5, col.lm = "green", lwd=par("lwd"), ...)
{
	#print(paste("Plotting xylm..."));
	#points(x, y, pch = 20, col = col, bg = bg, cex = cex);
	smoothScatter(x, y, nrpoints=0, bg=bg, cex=cex);
	
	ok <- is.finite(x) & is.finite(y)
	if (any(ok)) {
		abline(lm(y~x,subset=ok), col = col.lm, ...);
		#abline(lowess(y~x,subset=ok), col = col.lm, ...);
	}
}

# now use the function for the epi data. (see figure)

#pairs(epi, lower.panel=panel.smooth, upper.panel=panel.cor)


#use the same function for yet another graph. 
#pairs(bfi, lower.panel=panel.smooth, upper.panel=panel.cor)
