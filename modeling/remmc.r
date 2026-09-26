like <- function(trace1, trace2, p.m.o, p.m.n, p.n.o, p.n.n) {
    n.match <- sum(trace1 == trace2, na.rm=T)
    n.mis <- sum(trace1 != trace2, na.rm=T)
    
    return ((p.m.o / p.m.n)^n.match * (p.n.o / p.n.n)^n.mis)
}

N <- 40         # Number of studied items
n.f <- 30       # Number of features per item
u <- .6         # Probability of storing a feature in a given study event
u.x <- 0.2       # Probability of any stored feature being drawn from context (rather than content)
p.x <- u.x      # System's estimated probability that a given feature was stored from context
g <- .5         # Base-rate probability of any particular feature value
c <- .8         # Probability of correctly storing a particular feature value
d <- 1.0        # Probability that a context feature is shared between the study and test contexts

nsims <- 100    # Number of MC simulations to run

p.m.o <- p.x * (d * (c + (1 - c) * g) + (1 - d) * g) + (1 - p.x) * (c + (1 - c) * g)
p.m.n <- p.x * (d * (c + (1 - c) * g) + (1 - d) * g) + (1 - p.x) * g

p.n.o <- p.x * (d * (1 - g) * (1 - c) + (1 - d) * (1 - g)) + (1 - p.x) * (1 - c) * (1 - g)
p.n.n <- p.x * (d * (1 - g) * (1 - c) + (1 - d) * (1 - g)) + (1 - p.x) * (1 - g)

like.o <- c()
like.n <- c()
like.o.max <- c()
like.n.max <- c()

for (n in 1:nsims) {
    context <- rep(0, n.f)
    context[runif(n.f) < g] <- 1
    
    phys <- matrix(rep(0, N*n.f), nrow=N)
    phys[runif(N*n.f) < g] <- 1
    
    # Store study list traces
    traces <- matrix(rep(NA, N*n.f), nrow=N)
    for (t in 1:N) {
        for (f in 1:n.f) {
            if (runif(1) < u) {
                if (runif(1) < c) {
                    if (runif(1) < u.x) {
                        traces[t,f] <- context[f]
                    }
                    else {
                        traces[t,f] <- phys[t,f]
                    }
                }
                else {
                    if (runif(1) < g) {
                        traces[t,f] <- 1
                    }
                    else {
                        traces[t,f] <- 0
                    }
                }
            }
        }
    }
    
    # Context drift
    for (f in 1:n.f) {
        if (runif(1) > d) {
            if (runif(1) < g) {
                context[f] <- 1
            }
            else {
                context[f] <- 0
            }
        }
    }
    
    # Compute likelihood ratios
    for (t in 1:N) {
        temp.like <- c()
        probe <- phys[t,]
        choices <- runif(n.f)
        probe[choices < u.x] <- context[choices < u.x]
        for (j in 1:N) {temp.like <- c(temp.like, like(traces[j,], probe, p.m.o, p.m.n, p.n.o, p.n.n))}
        like.o <- c(like.o, mean(temp.like))
        like.o.max <- c(like.o.max, max(temp.like))
        
        temp.like <- c()
        probe <- rep(0, n.f)
        probe[runif(n.f) < g] <- 1
        choices <- runif(n.f)
        probe[choices < u.x] <- context[choices < u.x]
        for (j in 1:N) {temp.like <- c(temp.like, like(traces[j,], probe, p.m.o, p.m.n, p.n.o, p.n.n))}
        like.n <- c(like.n, mean(temp.like))
        like.n.max <- c(like.n.max, max(temp.like))
    }
}