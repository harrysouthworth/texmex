winterThresholded <- winter[winter$O3> median(winter$O3),]
nMissing <- dim(winter)[1] - dim(winterThresholded)[1]

# gpd fitting
padded <- as.data.frame(rbind(winterThresholded,
                              cbind(O3=runif(n=nMissing,min=0,max=min(winterThresholded$O3)),
                                    NO2=runif(n=nMissing,min=0,max=min(winterThresholded$NO2)),
                                    NO=runif(n=nMissing,min=0,max=min(winterThresholded$NO)),
                                    SO2=runif(n=nMissing,min=0,max=min(winterThresholded$SO2)),
                                    PM10=runif(n=nMissing,min=0,max=min(winterThresholded$PM10)))))


x.gpd <- migpd(padded,mqu=.7,  penalty="none")
r.gpd <- migpd(winter,mqu=.7,  penalty="none")
ref.gpd <- makeReferenceMarginalDistribution(x.gpd,r.gpd,whichNoChange=1)

depEstWithReferenceMargin <- mexDependence(x.gpd,which=1,dqu=0.7,referenceMargin=ref.gpd)
depEstWithOriginalMargin <-  mexDependence(r.gpd,which=1,dqu=0.7)

expect_equivalent(depEstWithReferenceMargin$dep$coefficients, depEstWithOriginalMargin$dep$coefficients)
expect_equivalent(depEstWithReferenceMargin$dep$Z, depEstWithOriginalMargin$dep$Z)
expect_equivalent(depEstWithReferenceMargin$dep$dth, depEstWithOriginalMargin$dep$dth)
expect_equivalent(depEstWithReferenceMargin$dep$dqu, depEstWithOriginalMargin$dep$dqu)
expect_equivalent(depEstWithReferenceMargin$dep$which, depEstWithOriginalMargin$dep$which)
expect_equivalent(depEstWithReferenceMargin$dep$conditioningVariable, depEstWithOriginalMargin$dep$conditioningVariable)
expect_equivalent(depEstWithReferenceMargin$dep$loglik, depEstWithOriginalMargin$dep$loglik)

expect_equivalent(depEstWithReferenceMargin$dep$margins[[1]],
                  depEstWithOriginalMargin$dep$margins[[1]])

expect_equivalent(depEstWithReferenceMargin$dep$constrain, depEstWithOriginalMargin$dep$constrain)
expect_equivalent(depEstWithReferenceMargin$dep$v, depEstWithOriginalMargin$dep$v)

set.seed(10)
predictionWithReferenceMargin <- predict(depEstWithReferenceMargin)
set.seed(10)
predictionwithOriginalMargins <- predict(depEstWithOriginalMargin)

expect_equivalent(predictionWithReferenceMargin,predictionWithReferenceMargin)
