# clt_simulation.R
# Run: in terminal -> Rscript clt_simulation.R
# Or open in VS Code with R extension and Source the file.

set.seed(2025) # reproducible

# Helper: overlay theoretical normal on a histogram
overlay_normal <- function(x, pop_mean, pop_sd, n, main="") {
  hist(x, breaks=40, probability=TRUE, main=main,
       xlab="Sample mean", ylim=c(0, max(density(x)$y, dnorm(pop_mean, pop_mean, pop_sd/sqrt(n)))) )
  # theoretical normal density
  curve(dnorm(x, mean=pop_mean, sd=pop_sd/sqrt(n)), add=TRUE, lwd=2)
}

# Main function: performs simulations and creates plots for one distribution
run_for_dist <- function(rfun, pop_params, pop_name, pop_mean, pop_sd, nrep=10000, pop_sample_size=1000, out_pdf=NULL) {
  # rfun: function(n, params) that returns n draws
  # pop_params: list passed to rfun using do.call
  # pop_mean, pop_sd: theoretical mean and sd of population
  if(!is.null(out_pdf)) pdf(out_pdf, width=11, height=8.5) # landscape letter-size
  par(mfrow=c(2,3), mar=c(4,4,2,1)) # grid for plots (6 slots)
  
  # 1) Population histogram (1000 draws)
  pop_samp <- do.call(rfun, c(list(n=pop_sample_size), pop_params))
  hist(pop_samp, breaks=40, probability=TRUE, main=paste(pop_name, "- population (n=1000)"),
       xlab="Value")
  curve(dnorm(x, mean=mean(pop_samp), sd=sd(pop_samp)), add=TRUE, lwd=2, col="darkgray")
  
  # Sample sizes to examine:
  ns <- c(5, 30, 50, 100)
  stored_means <- list()
  for(i in seq_along(ns)) {
    n <- ns[i]
    # Monte Carlo: compute 10000 sample means
    means <- replicate(nrep, mean(do.call(rfun, c(list(n=n), pop_params))))
    stored_means[[as.character(n)]] <- means
    
    # Plot histogram with theoretical normal overlay
    main_title <- paste(pop_name, "- n =", n)
    overlay_normal(means, pop_mean, pop_sd, n, main=main_title)
    # Print a small subtitle with empirical mean and sd
    m_emp <- mean(means); sd_emp<-sd(means)
    legend("topright", legend=c(sprintf("emp mean=%.3f", m_emp), sprintf("emp sd=%.3f", sd_emp)),
           bty="n")
  }
  # Leave last subplot empty if grid leftover
  par(mfrow=c(1,1))
  if(!is.null(out_pdf)) dev.off()
  
  # Return stored_means for further numerical checks if needed
  return(stored_means)
}

# Define distributions and their theoretical mean/sd
# 1) Uniform(0,1)
r_unif01 <- function(n) runif(n, 0, 1)
mu_unif <- 0.5; sd_unif <- sqrt(1/12)

# 2) Exponential(lambda = 1) -> rexp uses rate
r_exp1 <- function(n) rexp(n, rate = 1)
mu_exp <- 1; sd_exp <- 1

# 3) Binomial(n=10, p=0.5)
r_binom10 <- function(n) rbinom(n, size=10, prob=0.5)
mu_binom <- 10 * 0.5; sd_binom <- sqrt(10 * 0.5 * 0.5)

# 4) Chi-square(df=2) (highly skewed)
r_chisq2 <- function(n) rchisq(n, df = 2)
mu_chisq2 <- 2; sd_chisq2 <- sqrt(2*2)

# 5) t(df=3) (heavy-tailed)
r_t3 <- function(n) rt(n, df = 3)
# mean of t(df) is 0 for df>1, variance = df/(df-2) for df>2 -> here df=3: var=3/(3-2)=3 -> sd = sqrt(3)
mu_t3 <- 0
sd_t3 <- sqrt(3)

# 6) Normal(mean=4, sd=2)
r_norm42 <- function(n) rnorm(n, mean = 4, sd = 2)
mu_norm42 <- 4; sd_norm42 <- 2

# Create a PDF containing all distributions (one page per distribution)
out_pdf_all <- "CLT_simulations.pdf"
pdf(out_pdf_all, width=11, height=8.5)
dev.off() # create/overwrite empty file

# Run and save each distribution to its own page in the PDF (append)
run_and_append <- function(...) {
  tmp <- tempfile(fileext = ".pdf")
  r <- run_for_dist(..., out_pdf = tmp)
  # append tmp to main PDF (simple way: read pages copy -> but base R can't append easily)
  # Easier: we'll create separate files, then merge using system pdftk if available.
  # Instead, simply save each dist to its own file whose names you can print.
  # For convenience, we already saved each distribution's pdf in run_for_dist's out_pdf.
  invisible(r)
}

# We'll save separate PDFs per distribution so user can print them individually.
run_for_dist(r_unif01, list(), "Uniform(0,1)", mu_unif, sd_unif, out_pdf = "Uniform_0_1.pdf")
run_for_dist(r_exp1, list(), "Exponential(rate=1)", mu_exp, sd_exp, out_pdf = "Exponential_1.pdf")
run_for_dist(r_binom10, list(), "Binomial(n=10,p=0.5)", mu_binom, sd_binom, out_pdf = "Binomial_10_0.5.pdf")
run_for_dist(r_chisq2,_
