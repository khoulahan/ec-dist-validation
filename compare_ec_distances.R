### compare_ec_distances.R ########################################################################
# compare distance of evolutionary coupled residues in protein of interest versus queried
# population 
###################################################################################################
### PREAMBLE ######################################################################################
# load libraries
library(argparse);
library(BoutrosLab.plotting.general);
library(seqinr);

# general parameters
date <- Sys.Date();

### OBTAIN COMMAND LINE ARGUMENTS #################################################################
parser <- ArgumentParser();

parser$add_argument('-p', '--pdbIds', type = 'character', help = "filename of pdb IDs");
parser$add_argument('-t', '--threshold', type = 'integer', help = "EC score threshold");
parser$add_argument('-n', '--protein', type = 'integer', help = "ID for protein of interest");

args <- parser$parse_args();
###################################################################################################
# read in file of pdb Ids from query 
pdbIds <- read.delim(
    args$pdbIds,
    sep = '\t',
    header = FALSE
    );
# iterate over all ids and read in distance values 
metrics <- apply(
    pdbIds,
    1,
    function(id) {
        # read in distance file
        tmp <- read.delim(
            paste0(id, "_Distances.txt"),
            sep = '\t',
            header = FALSE
            );
        # extract only residues above specified EC score threshold
        # find median ec score, distance and number of residues above
        # the ec threshold 

        tmp <- tmp[tmp[,5] >= args$threshold,]
        # return metrics as dataframe
        data.frame(
            median = median(tmp[,6]),
            variance = var(tmp[,6]),
            group = 'PDB'
            );
metrics <- do.call('rbind', metrics);

### CALCULATE METRICS FOR PROTEIN OF INTEREST #####################################################
# read in distances for protein of interest
poi.distances <- read.delim(
    paste0(args$protein, "_Distances.txt"),
    sep = '\t',
    header = FALSE
    );
# only keep residues above specified EC score threshold
poi.distances   <- poi.distances[poi.distances[,5] > args$threshold,];
poi.metrics     <- data.frame(
    median = median(poi.distances[,6]),
    variance = var(poi.distances[,6]),
    group = 'POI'
    );

poi.percentile <- round(
    1-ecdf(metrics$median)(poi.metrics$median),
    digits = 2
    );

# combine poi with rest of metrics
metrics = rbind(metrics, poi.metrics)

### GENERATE SCATTERPLOT ##########################################################################
# order by median distance 
metrics <- metrics[order(metrics$median, decreasing = TRUE),];
metrics$index <-1:nrow(metrics);
# create scatterplot
create.scatterplot(
    formula = median ~ index,
    data = metrics,
    groups = metrics$group,
    filename = args$filename,
    yaxis.cex = 1.5,
    xaxis.cex = 1.5,
    ylab.label = 'Median Distance Between CA',
    ylab.cex = 1.5,
    xlab.label = 'Index',
    xlab.cex = 1.5,
    y.error.up = metrics$variance,
    error.bar.lwd = 1,
    error.whisker.angle = 120,
    y.error.bar.col = "black",
    col = c('black','firebrick3'),
    legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    points = list(
                        col = 'firebrick3',
                        pch = 21,
                        cex = 0.8,
                        fill = 'firebrick3'
                        ),
                    text = list(
                        lab = paste0(args$protein,': ', poi.percentile, ' Percentile')
                        ),
                    cex = 1
                    )
                ),
            x = 0.60,
            y = 0.97,
            draw = FALSE
            )
        )
    )