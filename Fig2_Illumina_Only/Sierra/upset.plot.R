library(UpSetR)

upset(fromList(listInput), order.by = "freq",mainbar.y.label="Number of peaks",nsets=7)
