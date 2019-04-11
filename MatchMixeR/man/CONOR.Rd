\name{MatchMixeR}
\alias{CONOR}
\alias{dwd}
\alias{eb}
\alias{xpn}
\alias{gq}

\title{
Functions for cross-platform normalization of microarray data, from the CONOR package.
}
\description{
The functions below perform cross-platform normalization of microarray data, and come from the CONOR package. Dwd is distance weighted discrimination. Eb is empirical bayes. Xpn is cross-platform normalization. The dwd function requires the package rdist.
The xpn function requires the package conclust.
}

\usage{

dwd(platform1.data, platform2.data, platform1.train = NULL,
    platform2.train = NULL, p1.names = 0, p2.names = 0, p1.train.names = 0,
    p2.train.names = 0, skip.match = FALSE, use.sparse = TRUE)

eb(platform1.data, platform2.data, par.prior = TRUE, filter = FALSE,
	prior.plots = FALSE, p1.names = 0, p2.names = 0, skip.match=FALSE)

xpn(platform1.data, platform2.data, K = 10, L = 4, p1.names = 0,
    p2.names = 0, gene.cluster = "kmeans", assay.cluster = "kmeans",
    corr = "pearson", iterations = 30, skip.match = FALSE)

gq(platform1.data, platform2.data, p1.names = 0, p2.names = 0, skip.match=FALSE)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{platform1.data}{
Expression data from platform 1.  Should be formatted as a data.frame, with each column representing an array and each row a gene.
}
  \item{platform2.data}{
Expression data from platform 2.  Should be formatted as a data.frame, with each column representing an array and each row a gene.
}
  \item{p1.names}{
One column of \code{platform1.data} may contain gene names.  The column number containing the gene names should be specified as \code{p1.names}.  If \code{p1.names} is zero, the rownames attribute will be used.  Default value is zero.  The gene names for platform 1 should correspond to the gene names for platform 2.
}
  \item{p2.names}{
One column of \code{platform2.data} may contain gene names.  The column number containing the gene names should be specified as \code{p2.names}.  If \code{p2.names} is zero, the rownames attribute will be used.  Default value is zero.  The gene names for platform 2 should correspond to the gene names for platform 1.
}

  \item{platform1.train}{Training data set for use with DWD.}
  \item{platform2.train}{Training data set for use with DWD.}
  \item{p1.train.names}{The column number containing the gene names for platform 1 training data.}
  \item{p2.train.names}{The column number containing the gene names for platform 2 training data.}
  \item{cluster}{Clustering method used by distran.}
  \item{par.prior}{Parameter for eb.  \code{par.prior} specifies whether to use a parametric or non-parametric prior distribution.}
  \item{filter}{Parameter for eb.}
  \item{prior.plots}{Generate prior plots for eb method.}


  \item{skip.match}{
If \code{skip.match} is \code{FALSE}, rows of \code{platform1.data} and \code{platform2.data} will be matched using gene names.  This process uses R's built-in \code{intercept} function, which can be quite inefficient.  If the rows of \code{platform1.data} and \code{platform2.data} already match, time can be saved by setting \code{skip.match} to \code{TRUE}.
}

  \item{use.sparse}{
For dwd only.  Can be set to \code{TRUE} or \code{FALSE}.  Determines whether dwd uses sparse matrices (via the Matrix package) for its internal calculations. Sparse matrix calculations are more efficient for large problems, but will not affect the final output.
}

  \item{gene.cluster}{
For \code{xpn} only, \code{gene.cluster} specifies the gene clustering method to be used.  Options are "kmeans", "pam", and "flexclust".  Only "kmeans" is practical for large numbers or genes.
}
  \item{assay.cluster}{
For \code{xpn} and \code{distran} only, \code{assay.cluster} specifies the assay clustering method to be used. Options are "classic", "kmeans", "pam", and "flexclust".
}
  \item{corr}{
For \code{xpn} and \code{distran} only, \code{corr} is the type of correlation to use as a distance measure for sample or gene clustering.  Ignored for the "kmeans" and "classic" clustering options, for which only Pearson's correlation is available.
}
  \item{iterations}{
For \code{xpn} only, \code{iterations} gives the number of iterations of the XPN algorithm to perform.
}

  \item{K}{
For \code{xpn} only, \code{K} is the number of gene clusters to use.  Must be an integer.  If a vector of more than one integer is given, the best value will be selected by the silhouette based method of the pamk function from the fpc package.  This can be extremely slow for large numbers of genes.
}

  \item{L}{
For \code{xpn} and \code{distran} only, \code{L} is the number of assay clusters to use.  Must be an integer.  If a vector of more than one integer is given, the best value will be selected by a silhouette based method using the pamk function from the fpc package.
}

}


\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{x}{Normalized data from platform 1.}
\item{y}{Normalized data from platform 2.}
\item{p1.adjust}{For dwd only, vector of platform effects removed from the platform 1 data.}
\item{p2.adjust}{For dwd only, vector of platform effects removed from the platform 2 data.}
}
\references{

Benito et al. Adjustment of systematic microarray data biases. Bioinformatics (2004) vol. 20 (1) pp. 105

Bolstad et al. A comparison of normalization methods for high density oligonucleotide array data based on variance and bias. Bioinformatics (2003) vol. 19 (2) pp. 185

Jiang et al. Joint analysis of two microarray gene-expression data sets to select lung adenocarcinoma marker genes. BMC bioinformatics (2004) vol. 5 pp. 81

Martinez et al. GenMiner: mining informative association rules from genomic data. Proceeding of the IEEE International Conference on Binformatics and Biomedicine. (2007) pp. 15-22

Shabalin et al. Merging two gene-expression studies via cross-platform normalization. Bioinformatics (2008) vol. 24 (9) pp. 1154

Shi et al. The MicroArray Quality Control (MAQC) project shows inter-and intraplatform reproducibility of gene expression measurements. Nature biotechnology (2006) vol. 24 (9) pp. 1151-1161

Walker et al. Empirical Bayes accomodation of batch-effects in microarray data using identical replicate reference samples: application to RNA expression profiling of blood .... BMC bioinformatics (2008)

Warnat et al. Cross-platform analysis of cancer microarray data improves gene expression based classification of phenotypes. BMC bioinformatics (2005) vol. 6 pp. 265

}
\author{
Jason Rudy and Faramarz Valafar
}



\examples{
fit_dwd <- dwd(gpl96,gpl570)
fit_eb <- eb(gpl96,gpl570)
fit_xpn <- xpn(gpl96,gpl570)
fit_gq <- gp(gpl96,gpl570)
}

\keyword{MatchMixeR}
\keyword{dwd}
\keyword{eb}
\keyword{xpn}
\keyword{distance weighted discrimination}
\keyword{empirical bayes}
\keyword{cross-platform normalization}

