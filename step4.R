# Authors: William John Fitzgerald & Tomas William Fitzgerald 
# main changepoint detection algorithm
# single changepoint model
step4 <- function(d) {
    n = length(d)
    dbar = mean(d)
    dsbar = mean(d*d)
    fac = dsbar-(dbar^2)
    summ = sum(d)
    summup = cumsum(d)
    y = vector()
    for (m in 1:(n-1)) {
        pos=m
        mscale = 4*(pos)*(n-pos)
        Q = summup[m]-(summ-summup[m])
        U = -(dbar*(n-2*pos) + Q)^2/mscale + fac
        y[m] = (-(n/2-1)*log(n*U/2) - 0.5*log((pos*(n-pos))))
    }
    return (c(y))
}

istep4 <- function(dat, np = 0.05) {
    rindexes = NULL
    indexes = c(1, length(dat))
    orders = NULL
    liklihoods = c(NA,NA)
    c = 0
    while(c<length(dat)) {
        chgs = vector()
        likli = vector()
        for(x in 1:(length(indexes)-1)) {
            if(length(dat[indexes[x]:(indexes[x+1])]) > 5) {
                steps = step4(dat[(indexes[x]+1):(indexes[x+1]-1)])
                likli[x] = max(steps)
                chgs[x] = which.max(steps)+indexes[x]+1
            }
        }
        indexes = c(indexes, chgs[complete.cases(chgs)])
        liklihoods = c(liklihoods,likli[complete.cases(chgs)])
        rindexes = c(rindexes,chgs[complete.cases(chgs)])
        liklihoods = liklihoods[order(indexes)]
        indexes = indexes[order(indexes)]
        orders = c(orders, rep(c, length(chgs[complete.cases(chgs)])))
        if(length(indexes)>length(dat)*np) {
            break
        }
        c = c+1
    }
    return(indexes)
}

trim <- function(indexes, t=10) {
    di =diff(indexes)
    while(min(di)<t) {
        indexes = indexes[-((which(di<t)[1]+1))]
        di =diff(indexes)
    }
    return(indexes)
}

seg_means <- function(data, indexes) {
    me = vector()
    for(x in 1:(length(indexes)-1)) {
        me[x] = mean(data[indexes[x]:(indexes[x+1])])
    }
return(me)
}

d = read.table("~/scratch/nanopore/10/A_0000018915.signal")[,1]
p = read.table("~/scratch/nanopore/10/A_0000018915.label")



indexes = istep4(d)
tindexes = trim(indexes, 10)
mmeans = seg_means(d, tindexes)
mdiffs = c(NA, diff(mmeans), NA)

plot(d)
for(x in 1:length(mmeans)) {
    segments(tindexes[x], mmeans[x], tindexes[x+1], mmeans[x], col="red")
}
