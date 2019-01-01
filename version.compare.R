# adapted from https://github.com/andrie/version.compare


library(version.compare)
library(knitr)
scale.factor <- 1.0   # scales the test sets down to x% to speed up the process for illustration

# r <- switch(Sys.info()[["sysname"]],
#             Linux = {
#               rscript <- findRscript()
#               
#               rv <- version.time(rscript, {
#                 as.character(getRversion())
#               })
#               idx <- which(unlist(rv$results) == "3.4.1")
#               rscript[idx]
#               
#             },
#             Windows = findRscript(version = "3.4.1.*x64"
#             )
# )


R.version
# Configure which installed version to use
r <- switch(
  Sys.info()[["sysname"]],
  Windows = c(
    "c:/program files/Microsoft/R Open/R-3.5.1/bin/x64/Rscript.exe", 
    "c:/program files/R/R-testing/bin/x64/Rscript.exe",
    "c:/program files/R/R-3.5.1/bin/x64/Rscript.exe"
  )
  # ,
  # Linux = c(
  #   # "/usr/lib64/RRO-8.0/R-3.1.1/lib/R/bin/Rscript",
  #   # "/usr/lib/R/bin/Rscript"
  # )
)


test.results <- RevoMultiBenchmark(rVersions = r, 
                                   threads = c(1, 8), # c(1, 4, 8), c(1, 8, 16),   # 
                                   scale.factor = scale.factor)

## C:\Program Files\Microsoft\R Open\R-3.4.1\bin\x64\Rscript.exe
## C:\Program Files\R\R-3.4.1\bin\x64\Rscript.exe

kable(test.results)
plot(test.results, theme_size = 8, main = "Elapsed time")

kable(urbanekPerformance(test.results), digits = 2)
plot(urbanekPerformance(test.results), theme_size = 8, main = "Relative Performance")



###################################################################################################

# Find installed versions of Rscript
#rscript <- findRscript(version = "3.4.1.*x64"


# # Configure which installed version to use
# r <- switch(
#   Sys.info()[["sysname"]],
#   Windows = c(
# #    "c:/program files/Microsoft/R Open/R-3.5.1/bin/x64/Rscript.exe", 
#     "c:/program files/R/R-testing/bin/x64/Rscript.exe",
#     "c:/program files/R/R-3.5.1/bin/x64/Rscript.exe"
#   )
#   # ,
#   # Linux = c(
#   #   # "/usr/lib64/RRO-8.0/R-3.1.1/lib/R/bin/Rscript",
#   #   # "/usr/lib/R/bin/Rscript"
#   # )
# )
# 
# 
# # Compute vector mean in different R installations
# version.time({
#   foo <- rnorm(1e9)   # rnorm(1e6)
#   mean(foo)
# } , r)
# 
# 
# # Compute matrix cross product in different R installations
# version.time({
#   m <- matrix(runif(100), nrow=10)
#   crossprod(m)
# } , r)

