# This code reproduces Figure 1.
rm(list = ls())
library(latex2exp)
library(energy)
library(ggplot2)
library(msda)
source("utility.R")
source("seas.R")
data("lymphoma", package = "spls")
#############################################
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

x <- lymphoma$x
y <- lymphoma$y
y <- y + 1

# Construct dummy variable for y.
y_dummy <- sapply(unique(y), function(i){
  as.numeric(y == i)
})

# Variable screen: keep 200 variables.
dist_cor <- sapply(seq_len(dim(x)[2]), function(i){
  dcor(y_dummy, x[,i])
})
ord <- order(dist_cor, decreasing = TRUE)[1:200]
x <- x[,ord, drop=FALSE]

######## Estimation consistency ##########
foo <- function(x, y){
  seassir_fit <- cv.seas(x, y, categorical=TRUE, nfolds = 5, lam1_fac=seq(1.2,0.5, length.out = 10), type = 'sir')
  beta_seassir <- seassir_fit$beta
  if(is.vector(beta_seassir)) beta_seassir <- matrix(beta_seassir)
  if(is.null(beta_seassir) || (seassir_fit$rank==0)){
    d_seassir <- NA
    ind_hat_seassir <- NA
    s_seassir <- NA
  }else{
    d_seassir <- seassir_fit$rank
    ind_hat_seassir <- which(apply(beta_seassir, 1, function(x) any(x!=0)))
    s_seassir <- length(ind_hat_seassir)
  }
  
  output <- list(rank = d_seassir, s = s_seassir, ind_hat = ind_hat_seassir, beta = beta_seassir)
  output
}

output <- foo(x,y) # One replicate
beta_seassir <- output$beta # SEAS-SIR
newx_seassir <- x %*% beta_seassir # Reduced predictor from SEAS-SIR

# Make plot using ggplot package.
data <- data.frame(x = newx_seassir[,1], y = newx_seassir[,2], class = factor(y))
g <- ggplot(data, aes(x=x, y=y))+
  geom_point(aes(shape = class), size = 3)+
  scale_shape_discrete(labels = c("DLBCL", "FL", "CLL"))+
  xlab(TeX('$\\hat{\\beta}_1^T\\mathbf{X}$'))+
  ylab(TeX('$\\hat{\\beta}_2^T\\mathbf{X}$'))+
  theme_bw()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 14),
        title = element_text(size = 14),
        legend.text = element_text(size = 16),
        legend.title = element_blank())
print(g)

## Colorful version.
# data <- data.frame(x = newx_seassir[,1], y = newx_seassir[,2], class = factor(y))
# g <- ggplot(data, aes(x=x, y=y))+
#   geom_point(aes(color = class), size = 3)+
#   scale_color_discrete(labels = c("DLBCL", "FL", "CLL"))+
#   xlab(TeX('$\\hat{\\beta}_1^T\\mathbf{X}$'))+
#   ylab(TeX('$\\hat{\\beta}_2^T\\mathbf{X}$'))+
#   theme_bw()+
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 14),
#         title = element_text(size = 14),
#         legend.text = element_text(size = 16),
#         legend.title = element_blank())
# print(g)

# Figure size: 600x400