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

# plot_heat_map(df, plot_title="Feature Table")
plot_heat_map <- function(df, plot_title="", pdf_outfile=character(0))
{
	data_matrix <- as.matrix(df);
	#data_matrix <- t(data_matrix);
	data_matrix <- t(data_matrix[nrow(data_matrix):1,]);
	
	cnames <- colnames(df);
	rnames <- rev(rownames(df));
	
	# Load RColorBrew
	library(RColorBrewer);
	pal <- brewer.pal(9, "Blues");
	

	
	
	if(length(pdf_outfile) != 0)
	{
		print(paste("Plotting ", plot_title, " to file ", pdf_outfile, sep=""));
		pdf_width <- ((ncol(df) * 40) + 200) / 75;
		pdf_height <- ((nrow(df) * 20) + 200) / 75;	
		pdf(pdf_outfile, width=pdf_width, height=pdf_height);
		
	} else {
		print(paste("Plotting ", plot_title, sep=""));
	}
	
	
	#
	# option	description
	# mar	numerical vector indicating margin size c(bottom, left, top, right) in lines. default = c(5, 4, 4, 2) + 0.1
	# mai	numerical vector indicating margin size c(bottom, left, top, right) in inches
	# pin	plot dimensions (width, height) in inches
	# http://research.stowers-institute.org/efg/R/Graphics/Basics/mar-oma/index.htm
	par(mar=c(3,14,15,2), oma=c(0.2,0.2,0.2,0.5), mex=0.5, mai=c(0.1,1, 2, 0.1));
	#par(mar=c(3,14,19,2), oma=c(0.2,0.2,0.2,0.5), mex=0.5);
	
	
	image(x=1:nrow(data_matrix), y=1:ncol(data_matrix), z=data_matrix, col=brewer.pal(6, "Blues"), ylab="", xlab="", axes=F);
	#image(x=1:nrow(data_matrix), y=(1+group_label_height):(ncol(data_matrix) + group_label_height), z=data_matrix, col=brewer.pal(6, "Blues"), ylab="", xlab="", axes=F);
	
	abline(h=c(1:ncol(data_matrix))+0.5 , v=c(1:nrow(data_matrix))+0.5, col="white", lwd=2, xpd=F);
	#abline(h=c((1:(ncol(data_matrix))+group_label_height))+0.5 , v=c(1:nrow(data_matrix))+0.5, col="white", lwd=2, xpd=F);
	
	
	# Y-axis labels
	axis(side=2, at=1:ncol(data_matrix), labels=rnames, col="white", las=1, cex.axis=0.85);
	#axis(side=2, at=((1+group_label_height):(ncol(data_matrix)+group_label_height)), labels=rnames, col="white", las=1, cex.axis=0.85);
	
	
	# X-axis labels
	text(1:nrow(data_matrix), par("usr")[4] + 0.3,  srt=55, adj=0, labels=cnames, xpd=T, cex=0.85);
	
	# Main Title
	text(par("usr")[1], par("usr")[4] + 5, plot_title, xpd=T, font=2, cex=1.5);
	
	if(length(pdf_outfile) != 0)
	{
		dev.off();
	} 
}






# 
test_code2 <- function() 
{

	df <- data.frame(matrix(sample(100, 600, replace=T), 40, 15));
	colnames(df) <- paste("Individual ", seq(1, ncol(df)), sep="");
	rownames(df) <- paste("Feature ", seq(1, nrow(df)), sep="");

	plot_heat_map(df, plot_title="Feature Table")



	plot_heat_map(df, plot_title="Feature Table", pdf_outfile="test.pdf");
	system("convert test.pdf test.png");


	df <- data.frame(matrix(sample(100, 100, replace=T), 20, 5));
	colnames(df) <- paste("Individual ", seq(1, ncol(df)), sep="");
	rownames(df) <- paste("Feature ", seq(1, nrow(df)), sep="");
}



# grid.arrange draws directly on a device. arrangeGrob, on the other hand, doesn't draw anything but returns a grob g, that you can pass to ggsave(file="whatever.pdf", g).		
# The reason it works differently than with ggplot objects, where by default the last plot is being saved if not specified, is that ggplot2 keeps a track of the last plot,and I don't think grid.arrange should mess with this counter private to the package.
# http://stackoverflow.com/questions/18427455/multiple-ggplots-of-different-sizes
test_code_heatmap_by_group <- function()
{
	source("/data/bcnskaa/1KGP/Small_Indels/process_DF.R");
	source("/data/bcnskaa/1KGP/Small_Indels/plot_heat_map.R");
	
	SELECTED_COLOR_SET = c("Red", "Orangered", "Orange", "Yellow", "Yellow3", "Green", "Olivedrab", "Darkgreen", "Darkcyan", "Blue", "Darkslateblue", "Darkviolet");
	SELECTED_COLOR_SET2 = c("lightpink", "paleturquoise", "moccasin", "lightblue1", "palegreen");

	row_n <- 20;
	col_n <- 15;
	
	x_labels <- paste("Feature_", seq(1, col_n), sep="");
	x_grouplabels <- sample(c("Group1", "Group2", "Group3"), length(x_labels), replace=T);
	y_labels <- paste("Individual_", letters[1 : row_n], sep="");
	y_grouplabels <- sample(c("A", "B", "C"), length(y_labels), replace=T);

	
	xlabel_meta_info <- data.frame(LABEL=x_labels, GROUP=x_grouplabels, stringsAsFactors=F);
	ylabel_meta_info <- data.frame(LABEL=y_labels, GROUP=y_grouplabels, stringsAsFactors=F);
	
	
	z <- sample(seq(0,1, by=0.01), row_n * col_n, replace=T)
	
	
	plot_mtx <- data.frame(matrix(z, length(y_labels), length(x_labels)), stringsAsFactors=F);
	colnames(plot_mtx) <- x_labels;
	rownames(plot_mtx) <- y_labels;
	
	
	plot_data <- create_plot_data_from_plot_mtx(plot_mtx);
	source("/data/bcnskaa/1KGP/Small_Indels/plot_heat_map.R")

	
	plot_heatmap_with_group_info(plot_data, xlabel_meta_info, ylabel_meta_info, plot_title="Test Correlation");
	
}


# http://stackoverflow.com/questions/13294952/left-align-two-graph-edges-ggplot
# label_meta_info: [LABEL, GROUP, COLOR]
# 
plot_heatmap_with_group_info <- function(plot_data, xlabel_meta_info=data.frame(), ylabel_meta_info=data.frame(), plot_title=character(0), show_value=F, pdf_fn=character(0), pdf_width=8, pdf_height=8, title_x=-5, title_y=1.05, scale_limits=character(0), scale_colors=character(0)) 
{
	SELECTED_COLOR_SET = c("Red", "Orangered", "Orange", "Yellow", "Yellow3", "Green", "Olivedrab", "Darkgreen", "Darkcyan", "Blue", "Darkslateblue", "Darkviolet");
	SELECTED_COLOR_SET2 = c("olivedrab1", "paleturquoise", "moccasin", "gold3", "brown2", "mediumorchid", "lightpink", "lightsalmon");
	
	library(ggplot2);
	library(RColorBrewer);
	library(grid);
	library(gridExtra);
	library(scales);
	source("/data/bcnskaa/1KGP/Small_Indels/process_DF.R");

	
	if(length(scale_limits) == 0)
	{
		scale_limits <- c(-1.0,1.0)	
	}

	if(length(scale_colors) == 0)
	{
		scale_colors <- c("red3", "white", "steelblue");
	}
	
	# A plot
	#labels <- paste("Feature_", seq(1,25), sep="")
	#groups <- paste("Group", sample(c(1,2,3,4), length(labels), replace=T), sep="")
	#x <- rep(1, length(labels));
	#y <- seq(1, length(labels));
	#plot_data <- data.frame(x=x, y=y, label=labels, group=groups, stringsAsFactors=F);
	
	
	xseparators <- data.frame(X=integer(0), Y=integer(0));
	xgrouplabels <- data.frame(X=integer(0), Y=integer(0), GROUP=character(0));
	yseparators <- data.frame(X=integer(0), Y=integer(0));
	ygrouplabels <- data.frame(X=integer(0), Y=integer(0), GROUP=character(0));
	
	
	# The basic rule with ggplot2 that applies to almost anything that you want in a specific order is: If you want something to appear in a particular order, you must make the corresponding variable a factor, with the levels sorted in your desired order.
	# Prepare group label for x-axis labels
	if(dim(xlabel_meta_info)[1] != 0)
	{
		plot_data$X <- factor(plot_data$X, levels=xlabel_meta_info$LABEL, ordered=T);
		
		labels <- xlabel_meta_info$LABEL
		groups <- xlabel_meta_info$GROUP
		
		# Generate separators for group labels
		xseparators = data.frame(X=integer(0), Y=integer(0))
		for(i in 1 : (length(groups) - 1))
		{
			g1 <- groups[i]
			g2 <- groups[i+1]
			
			if(g1 != g2)
			{
				# add a separator
				xseparators <- rbind(xseparators, data.frame(Y=-0.5, X=i+0.5))
			}
		}

		xgrouplabels <- data.frame(Y=0, X=xlabel_meta_info$LABEL, GROUP=xlabel_meta_info$GROUP);
		xgrouplabels$X <- factor(xgrouplabels$X, levels=xlabel_meta_info$LABEL, ordered=T)
		
		# Colors for x label groups
		x_label_group_colors <- SELECTED_COLOR_SET2[1:length(unique(as.character(xgrouplabels$GROUP)))]
	}
	
	
	# Prepare group label for y-axis labels
	if(dim(ylabel_meta_info)[1] != 0)
	{
		plot_data$Y <- factor(plot_data$Y, levels=ylabel_meta_info$LABEL, ordered=T);
		
		labels <- ylabel_meta_info$LABEL
		groups <- ylabel_meta_info$GROUP
		
		# Generate separators for group labels
		yseparators = data.frame(X=integer(0), Y=integer(0))
		for(i in 1 : (length(groups) - 1))
		{
			g1 <- groups[i]
			g2 <- groups[i+1]
			
			if(g1 != g2)
			{
				# add a separator
				yseparators <- rbind(yseparators, data.frame(X=-0.5, Y=i+0.4))
			}
		}
		
		ygrouplabels <- data.frame(X=0, Y=ylabel_meta_info$LABEL, GROUP=ylabel_meta_info$GROUP);
		ygrouplabels$Y <- factor(ygrouplabels$Y, levels=ylabel_meta_info$LABEL, ordered=T)
		
		# Colors for y label groups
		y_label_group_colors <- SELECTED_COLOR_SET2[1:length(unique(as.character(ygrouplabels$GROUP)))]
	}

	
	
	# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

	# y group label
	yg <- ggplot(data=ygrouplabels) + geom_tile(data=ygrouplabels, aes(x=X, y=Y, fill=GROUP, height=1.2, width=0.3)) + 
		scale_fill_manual(values=y_label_group_colors) +
		geom_segment(data=yseparators, aes(xend=X+2, x=-1, y=Y, yend=Y), size=1, color="white") +
		theme(axis.text.x=element_blank(), axis.text.y=element_text(hjust=1, vjust=0.5)) +
		theme(plot.margin = unit(c(0,-0.3,0,0), "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank());

	ygroup_legend <- g_object(yg);

	
	# x group label
	xg <- ggplot(data=xgrouplabels) + geom_tile(data=xgrouplabels, aes(x=X, y=Y, fill=GROUP, height=0.3, width=1)) + 
			scale_fill_manual(values=x_label_group_colors) +
			geom_segment(data=xseparators, aes(yend=Y+2, y=-1, x=X, xend=X), size=1, color="white") +
			#theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
			theme(axis.text.y=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
			theme(plot.margin = unit(c(-2,0, 1,-1), "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank());

	xgroup_legend <- g_object(xg);
	
	

	# Main 
	g_main <- ggplot() + geom_tile(aes(x=X, y=Y, fill=Z, height=0.97, width=0.97), data=plot_data) +
			#scale_fill_gradient(low="Red3", high="SteelBlue") + 
			#scale_fill_gradient2(low="red3", high="steelblue") +
			scale_fill_gradientn(colours=scale_colors, limits=scale_limits) +
			#geom_text(aes(x=X,y=Y,label=Z), size=3, data=plot_data) + 
			#scale_y_discrete(expand=c(0,1,1,1)) +
			scale_y_discrete() + 
			#geom_text(aes(x=nr_xlabels,y=rep(7.2, length(nr_xlabels)), label=nr_xlabels, angle=75),size=4) +
			
			# http://stackoverflow.com/questions/7263849/what-do-hjust-and-vjust-do-when-making-a-plot-using-ggplot
			theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
			# theme(plot.margin=unit(c(1,1,1,1), "cm")) +
			#theme(axis.text.x=element_blank()) + 
			theme(plot.margin = unit(c(0,0,0,-1.0), "cm"), panel.background=element_blank(),  axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
	
	
	if(show_value == T)
	{
		g_main <- g_main + geom_text(aes(x=X,y=Y,label=Z), size=3, data=plot_data);
	}
	


	# top, right, bottom, and left
	#theme(plot.margin=unit(c(0.5,0.1,0.1,0.1), "cm"));
	
	# legend.position="none",
	main_legend <- g_object(g_main)
	
	# Final touch ups
	g_main <- g_main + labs(x=NULL, y=NULL) + theme(axis.text=element_blank(), legend.position="none");
	yg <- yg + theme(legend.position = 'none') + labs(x=NULL);
	xg <- xg + theme(legend.position = 'none') + labs(y=NULL);
	
	
	
	# http://stackoverflow.com/questions/12539348/ggplot-separate-legend-and-plot

	print("Plotting graph...");
	if(length(pdf_fn) != 0)
	{
		pdf(pdf_fn, width=pdf_width,height=pdf_height);
	}
	
	
	# http://stackoverflow.com/questions/8112208/how-can-i-obtain-an-unbalanced-grid-of-ggplots
	# http://stackoverflow.com/questions/18427455/multiple-ggplots-of-different-sizes
	# http://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
	vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)	


	grid.newpage();
	pushViewport(viewport(layout = grid.layout(12, 12))) # 10 rows, 11 columns

	#print(ggplot() + ggtitle("Hello"), vp=vplayout(1,1:12))
	print(yg, vp=vplayout(2:10,1:2))
	print(xg, vp=vplayout(11:12,3:11))
	print(g_main, vp=vplayout(2:10,3:11))
	
	
	pushViewport(vp=vplayout(2:12,12))
	g_legends <- arrangeGrob(ygroup_legend, main_legend, widths=c(1/2,1/2), nrow=2)
	grid.draw(g_legends)
	
	if(length(plot_title) != 0)
	{
		grid.text(plot_title, x=title_x, y=title_y, vp=vplayout(2:10,1:12)); 
		#grid.text(plot_title, x=-5, y=1.05, vp=vplayout(1,1:12)); 
	}
	

	

#	# http://stackoverflow.com/questions/8112208/how-can-i-obtain-an-unbalanced-grid-of-ggplots
#	# http://stackoverflow.com/questions/18427455/multiple-ggplots-of-different-sizes
#	vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#	grid.newpage();
#	pushViewport(viewport(layout = grid.layout(12, 12))) # 10 rows, 11 columns
#	
#
#	#print(ggplot() + ggtitle("Hello"), vp=vplayout(1,1:12))
#	print(yg, vp=vplayout(2:10,1:2))
#	print(xg, vp=vplayout(11:12,3:11))
#	print(g_main, vp=vplayout(2:10,3:11))
#	#g_legends <- arrangeGrob(ygroup_legend, legend, widths=c(1/2,1/2), nrow=2)
#	print(legend, vp=vplayout(1:12,12))
#
#	
	


	if(length(pdf_fn) != 0)
	{
		dev.off();

	}

}



# http://stackoverflow.com/questions/13294952/left-align-two-graph-edges-ggplot
# label_meta_info: [LABEL, GROUP, COLOR]
# 
plot_heatmap_with_double_y <- function(plot_data, xylabel_meta_info=data.frame(), ylabel_meta_info=data.frame(), plot_title=character(0), show_value=F, pdf_fn=character(0), pdf_width=8, pdf_height=8, title_x=-5, title_y=1.05, scale_limits=character(0), scale_colors=character(0)) 
{
	SELECTED_COLOR_SET = c("Red", "Orangered", "Orange", "Yellow", "Yellow3", "Green", "Olivedrab", "Darkgreen", "Darkcyan", "Blue", "Darkslateblue", "Darkviolet");
	SELECTED_COLOR_SET2 = c("olivedrab1", "paleturquoise", "moccasin", "gold3", "brown2", "mediumorchid", "lightpink", "lightsalmon");
	
	library(ggplot2);
	library(RColorBrewer);
	library(grid);
	library(gridExtra);
	library(scales);
	source("/data/bcnskaa/1KGP/Small_Indels/process_DF.R");
	
	
	if(length(scale_limits) == 0)
	{
		scale_limits <- c(-1.0,1.0)	
	}
	
	if(length(scale_colors) == 0)
	{
		scale_colors <- c("red3", "white", "steelblue");
	}
	

	xseparators <- data.frame(X=integer(0), Y=integer(0));
	xgrouplabels <- data.frame(X=integer(0), Y=integer(0), GROUP=character(0));
	yseparators <- data.frame(X=integer(0), Y=integer(0));
	ygrouplabels <- data.frame(X=integer(0), Y=integer(0), GROUP=character(0));
	
	
	# The basic rule with ggplot2 that applies to almost anything that you want in a specific order is: If you want something to appear in a particular order, you must make the corresponding variable a factor, with the levels sorted in your desired order.
	# Prepare group label for x-axis labels
	if(dim(xlabel_meta_info)[1] != 0)
	{
		plot_data$X <- factor(plot_data$X, levels=xlabel_meta_info$LABEL, ordered=T);
		
		labels <- xlabel_meta_info$LABEL
		groups <- xlabel_meta_info$GROUP
		
		# Generate separators for group labels
		xseparators = data.frame(X=integer(0), Y=integer(0))
		for(i in 1 : (length(groups) - 1))
		{
			g1 <- groups[i]
			g2 <- groups[i+1]
			
			if(g1 != g2)
			{
				# add a separator
				xseparators <- rbind(xseparators, data.frame(Y=-0.5, X=i+0.5))
			}
		}
		
		xgrouplabels <- data.frame(Y=0, X=xlabel_meta_info$LABEL, GROUP=xlabel_meta_info$GROUP);
		xgrouplabels$X <- factor(xgrouplabels$X, levels=xlabel_meta_info$LABEL, ordered=T)
		
		# Colors for x label groups
		x_label_group_colors <- SELECTED_COLOR_SET2[1:length(unique(as.character(xgrouplabels$GROUP)))]
	}
	
	
	# Prepare group label for y-axis labels
	if(dim(ylabel_meta_info)[1] != 0)
	{
		plot_data$Y <- factor(plot_data$Y, levels=ylabel_meta_info$LABEL, ordered=T);
		
		labels <- ylabel_meta_info$LABEL
		groups <- ylabel_meta_info$GROUP
		
		# Generate separators for group labels
		yseparators = data.frame(X=integer(0), Y=integer(0))
		for(i in 1 : (length(groups) - 1))
		{
			g1 <- groups[i]
			g2 <- groups[i+1]
			
			if(g1 != g2)
			{
				# add a separator
				yseparators <- rbind(yseparators, data.frame(X=-0.5, Y=i+0.4))
			}
		}
		
		ygrouplabels <- data.frame(X=0, Y=ylabel_meta_info$LABEL, GROUP=ylabel_meta_info$GROUP);
		ygrouplabels$Y <- factor(ygrouplabels$Y, levels=ylabel_meta_info$LABEL, ordered=T)
		
		# Colors for y label groups
		y_label_group_colors <- SELECTED_COLOR_SET2[1:length(unique(as.character(ygrouplabels$GROUP)))]
	}
	

	
	# y group label
	yg <- ggplot(data=ygrouplabels) + geom_tile(data=ygrouplabels, aes(x=X, y=Y, fill=GROUP, height=1.2, width=0.3)) + 
			scale_fill_manual(values=y_label_group_colors) +
			geom_segment(data=yseparators, aes(xend=X+2, x=-1, y=Y, yend=Y), size=1, color="white") +
			theme(axis.text.x=element_blank(), axis.text.y=element_text(hjust=1, vjust=0.5)) +
			theme(plot.margin = unit(c(0,-0.3,0,0), "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank());
	
	ygroup_legend <- g_object(yg);
	
	
	# x group label
	xg <- ggplot(data=xgrouplabels) + geom_tile(data=xgrouplabels, aes(x=X, y=Y, fill=GROUP, height=0.3, width=1)) + 
			scale_fill_manual(values=x_label_group_colors) +
			geom_segment(data=xseparators, aes(yend=Y+2, y=-1, x=X, xend=X), size=1, color="white") +
			#theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
			theme(axis.text.y=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
			theme(plot.margin = unit(c(-2,0, 1,-1), "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank());
	
	xgroup_legend <- g_object(xg);
	
	
	
	# Main 
	g_main <- ggplot() + geom_tile(aes(x=X, y=Y, fill=Z, height=0.97, width=0.97), data=plot_data) +
			#scale_fill_gradient(low="Red3", high="SteelBlue") + 
			#scale_fill_gradient2(low="red3", high="steelblue") +
			scale_fill_gradientn(colours=scale_colors, limits=scale_limits) +
			#geom_text(aes(x=X,y=Y,label=Z), size=3, data=plot_data) + 
			#scale_y_discrete(expand=c(0,1,1,1)) +
			scale_y_discrete() + 
			#geom_text(aes(x=nr_xlabels,y=rep(7.2, length(nr_xlabels)), label=nr_xlabels, angle=75),size=4) +
			
			# http://stackoverflow.com/questions/7263849/what-do-hjust-and-vjust-do-when-making-a-plot-using-ggplot
			theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
			# theme(plot.margin=unit(c(1,1,1,1), "cm")) +
			#theme(axis.text.x=element_blank()) + 
			theme(plot.margin = unit(c(0,0,0,-1.0), "cm"), panel.background=element_blank(),  axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
	
	
	if(show_value == T)
	{
		g_main <- g_main + geom_text(aes(x=X,y=Y,label=Z), size=3, data=plot_data);
	}
	
	

	
	# top, right, bottom, and left
	#theme(plot.margin=unit(c(0.5,0.1,0.1,0.1), "cm"));
	
	# legend.position="none",
	main_legend <- g_object(g_main)
	
	# Final touch ups
	g_main <- g_main + labs(x=NULL, y=NULL) + theme(axis.text=element_blank(), legend.position="none");
	yg <- yg + theme(legend.position = 'none') + labs(x=NULL);
	xg <- xg + theme(legend.position = 'none') + labs(y=NULL);
	
	
	
	# http://stackoverflow.com/questions/12539348/ggplot-separate-legend-and-plot
	#grid.arrange(g + theme(legend.position = 'none') + labs(x=NULL), g_main + labs(x=NULL, y=NULL) + theme(axis.text=element_blank(), legend.position="none"), legend, ncol=3, widths=c(1/11, 9/11, 1/11))
	#grid.arrange(yg + theme(legend.position = 'none') + labs(x=NULL), g_main + labs(x=NULL, y=NULL) + theme(axis.text=element_blank(), legend.position="none"), legend, ncol=3, widths=c(1/11, 9/11, 1/11))
	#grid.arrange(yg + theme(legend.position = 'none') + labs(x=NULL), g_main + labs(x=NULL, y=NULL) + theme(axis.text=element_blank(), legend.position="none"), legend, nrow=3, ncol=3, widths=c(1/11, 9/11, 1/11))
	
	
	
	print("Plotting graph...");
	if(length(pdf_fn) != 0)
	{
		pdf(pdf_fn, width=pdf_width,height=pdf_height);
	}
	
	
	# http://stackoverflow.com/questions/8112208/how-can-i-obtain-an-unbalanced-grid-of-ggplots
	# http://stackoverflow.com/questions/18427455/multiple-ggplots-of-different-sizes
	# http://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
	vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)	
	
	
	grid.newpage();
	pushViewport(viewport(layout = grid.layout(12, 12))) # 10 rows, 11 columns
	
	#print(ggplot() + ggtitle("Hello"), vp=vplayout(1,1:12))
	print(yg, vp=vplayout(2:10,1:2))
	print(xg, vp=vplayout(11:12,3:11))
	print(g_main, vp=vplayout(2:10,3:11))
	
	
	pushViewport(vp=vplayout(2:12,12))
	g_legends <- arrangeGrob(ygroup_legend, main_legend, widths=c(1/2,1/2), nrow=2)
	grid.draw(g_legends)
	
	if(length(plot_title) != 0)
	{
		grid.text(plot_title, x=title_x, y=title_y, vp=vplayout(2:10,1:12)); 
		#grid.text(plot_title, x=-5, y=1.05, vp=vplayout(1,1:12)); 
	}
	

	
	if(length(pdf_fn) != 0)
	{
		dev.off();
		
	}
	
}





# http://stackoverflow.com/questions/13294952/left-align-two-graph-edges-ggplot
# label_meta_info: [LABEL, GROUP, COLOR]
# 
#plot_heatmap_group_info2 <- function(plot_data, xlabel_meta_info=data.frame(), ylabel_meta_info=data.frame(), plot_title=character(0), show_value=F, pdf_fn=character(0), pdf_width=8, pdf_height=8, scale_limits=character(0)) 
plot_heatmap <- function(plot_data, xlabel_meta_info=data.frame(), ylabel_meta_info=data.frame(), plot_title=character(0), pdf_fn=character(0), pdf_width=8, pdf_height=8, scale_limits=character(0), show_value=F) 
{
	SELECTED_COLOR_SET = c("Red", "Orangered", "Orange", "Yellow", "Yellow3", "Green", "Olivedrab", "Darkgreen", "Darkcyan", "Blue", "Darkslateblue", "Darkviolet");
	SELECTED_COLOR_SET2 = c("lightpink", "paleturquoise", "moccasin", "lightblue1", "palegreen");
	
	library(ggplot2);
	library(RColorBrewer);
	library(grid);
	library(scales);
	source("/data/bcnskaa/1KGP/Small_Indels/process_DF.R");
	
	scale_limits <- c(-1.0,1.0)

	
	if(dim(xlabel_meta_info)[1] != 0)
	{
		plot_data$X <- factor(plot_data$X, levels=as.character(xlabel_meta_info$LABEL), ordered=T)
	}
	
	if(dim(ylabel_meta_info)[1] != 0)
	{
		plot_data$Y <- factor(plot_data$Y, levels=as.character(ylabel_meta_info$LABEL), ordered=T)
	}
	

	if(length(pdf_fn) != 0)
	{
		pdf(pdf_fn, width=pdf_width,height=pdf_height);
	}
	
	# Main 
	g_main <- ggplot() + geom_tile(aes(x=X, y=Y, fill=Z, height=0.97, width=0.97), data=plot_data) +
			#scale_fill_gradient(low="Red3", high="SteelBlue") + 
			#scale_fill_gradient2(low="red3", high="steelblue") +
			scale_fill_gradientn(colours=c("red3", "white", "steelblue"), limits=scale_limits) +
			#geom_text(aes(x=X,y=Y,label=Z), size=3, data=plot_data) + 
			scale_y_discrete(expand=c(0,1,1,1)) + 
			#geom_text(aes(x=nr_xlabels,y=rep(7.2, length(nr_xlabels)), label=nr_xlabels, angle=75),size=4) +
			
			# http://stackoverflow.com/questions/7263849/what-do-hjust-and-vjust-do-when-making-a-plot-using-ggplot
			theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
			theme(plot.margin=unit(c(1,1,1,1), "cm")) +
			#theme(axis.text.x=element_blank()) + 
			theme(panel.background=element_blank(),  axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
	
	
	if(show_value == T)
	{
		g_main <- g_main + geom_text(aes(x=X,y=Y,label=Z), size=3, data=plot_data);
	}
	
	if(length(plot_title) != 0)
	{
		#g_main <- g_main + ggtitle(plot_title)
		g_main <- g_main + ggtitle(plot_title); 
	}
	
	print("Plotting graph...");
	print(g_main)
	
	if(length(pdf_fn) != 0)
	{
		dev.off();
		print(paste("Converting ", pdf_fn, "...", sep=""));
		system(paste("convert", pdf_fn, sub(".pdf", ".png", pdf_fn)));
	}
	
	# http://stackoverflow.com/questions/12539348/ggplot-separate-legend-and-plot
	#grid.arrange(g + theme(legend.position = 'none'), g_main + theme(axis.text=element_blank(), legend.position="none"), legend, ncol=3, nrow=2, widths=c(2/11, 7/11, 2/11))
}


test <- function()
{	
}


# Get the legend item from a ggplot object
g_object <- function(a.gplot, selected_name="guide-box")
{
	library(gridExtra)
	tmp <- ggplot_gtable(ggplot_build(a.gplot))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == selected_name)
	legend <- tmp$grobs[[leg]]
	legend
}





test_color_func <- function(z) {
	print(z)
	group_ids <- c("Group 1","Group 2","Group 3","Group 4","Group 5");
	col <- c("darkblue", "blue", "green", "orange", "red");
	return(unlist(col[match(z, group_ids)]));
}



plot_RColor_map <- function(pdf_fn=character(0)) {
	
# Print R internal colors
	available_colors <- sort(colors());
	available_colors <- available_colors[-grep("gray", available_colors)]
	
	color_n <- length(available_colors)
	col_n <- 20
	row_n <- ceiling(length(available_colors) / col_n)
	
	available_colors <- c(available_colors, rep("white", (col_n - (length(available_colors) %% col_n))))
	
	
	x <- rep(seq(1, col_n), row_n);
	y <- rep(seq(1, row_n),each=col_n)
	m <- matrix(seq(1,length(available_colors)), ncol=col_n)
	
	plot_mtx <- data.frame(m)
	plot_mtx[plot_mtx>color_n] <- 1
	rownames(plot_mtx) <- seq(1, row_n);
	colnames(plot_mtx) <- seq(1, col_n);
	
	
	source("/data/bcnskaa/1KGP/Small_Indels/process_DF.R")
	plot_data <- create_plot_data_from_plot_mtx(plot_mtx)
	plot_data <- cbind(plot_data, COLOR=available_colors)
	plot_data$X <- factor(plot_data$X, levels=seq(1, col_n), ordered=T);
	plot_data$Y <- factor(plot_data$Y, levels=seq(row_n, 1), ordered=T);
	
	
#guide_colors <- data.frame(X=plot_data$X,Y=plot_data$Y,ID=seq(1, color_n),LABEL=available_colors, COLOR=available_colors, stringsAsFactors=F)
	
	
	library(ggplot2)
	library(grid)
	
	
	if(length(pdf_fn) == 1)
	{
		pdf_fn = "R_color_map.pdf"
		pdf(pdf_fn, width=(col_n*1)+5, height=(row_n*1)+5)
	}
	
	ggplot() + geom_tile(aes(x=X, y=Y, fill=COLOR, height=0.97, width=0.97), data=plot_data) +
			ggtitle("R Color Codes") +
			scale_fill_manual(values=as.character(plot_data$COLOR)) +
			geom_text(aes(x=X,y=Y,label=COLOR), size=3, data=plot_data) + 
			scale_y_discrete(expand=c(0,1,1,1)) + 
			theme(axis.text.x=element_text(angle=0, hjust=0, vjust=0)) + 
			theme(plot.margin=unit(c(1,1,1,1), "cm")) +
			theme(legend.position="none", panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
	
	if(length(pdf_fn) == 1)
	{
		dev.off();
	}	
	
}
