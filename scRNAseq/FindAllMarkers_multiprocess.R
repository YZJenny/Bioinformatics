library(future)
# check the current active plan
plan()
plan("multiprocess", workers = 4)
plan()
start = Sys.time()
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
