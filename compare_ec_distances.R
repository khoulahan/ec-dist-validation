### compare_ec_distances.R ########################################################################
# compare distance of evolutionary coupled residues in protein of interest versus queried
# population 
###################################################################################################
### PREAMBLE ######################################################################################
# load libraries
library(argparse);

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
query.distances <- apply(
    pdbIds,
    1,
    function(id) {
        tmp <- read.delim(
            paste0(id, "_Distances.txt"),
            sep = '\t',
            header = FALSE
            );
        # extract only residues above specified EC score threshold
        # return only sum of distances that pass threshold 
        sum(tmp[tmp[,5] >= args$threshold,6])
        })
query.distances <- unlist(query.distances);

# read in distances for protein of interest
protein.distances <- read.delim(
    paste0(args$protein, "_Distances.txt"),
    sep = '\t',
    header = FALSE
    );
# only keep residues above specified EC score threshold
protein.distances   <- protein.distances[protein.distances[,5] > args$threshold,6];
protein.sum         <- sum(protein.distances);

# check that protein sum falls within range of pdb database
# if not, do not produce plot
if (protein.sum > max(query.distances)) {
    cat("Distances calculated for protein of interest larger than PDB distribution.\n");
    cat("Failed to produce density plot ...\n")
} else {
    # find protein percentile and round to nearest hundredth 
    protein.percentile <- round(
        1-ecdf(query.distances)(protein.sum),
        digits = 2
        );
    # generate density plot
    create.densityplot(
        x = data.frame(query.distances),
        filename = paste0(
            args$protein,
            '_DistancesRelativeToPDB.tiff'
            ),
        xlab.label = 'Sum of CA Distances',
        abline.v = protein.sum,
        abline.col = 'darkorchid4',
        legend = list(
            inside = list(
                fun = draw.key,
                args = list(
                    key = list(
                        points = list(
                            col = 'darkorchid4',
                            pch = 21,
                            cex = 0.8,
                            fill = 'darkorchid4'
                            ),
                        text = list(
                            lab = paste0(args$protein,': ', protein.percentile, " Percentile\n",length(protein.distances), " Residues Used")
                            ),
                        cex = 1
                        )
                    ),
                x = 0.50,
                y = 0.97,
                draw = FALSE
                )
            )
        );
    }