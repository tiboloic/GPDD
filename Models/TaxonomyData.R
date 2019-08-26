# LT 
# generate data structure and files for TMB fitting
# save in Rdata file

dataset = read.csv("GPDD2.csv")

list_length = list()     
list_series = list()
list_indices = list()

IDs = unique(dataset$ID)

# preraring the DAT file
for (ID in IDs) {
  
  iID = match(ID, IDs)
  
  
  # extract time series
  rawseries = dataset[dataset[,"ID"]==ID,]
  
  # problem with series GPPDD ID 1505 (iID=105)
  # only only obs every 5 steps
  # change step
  if (ID %in% c(1505, 1507, 1508, 1516, 1517)) { rawseries$SeriesStep = rawseries$SeriesStep %/% 5 }
  
  # replace 0 by NA. That's the best we can do with gaussian filter 
  # instead Poisson error should be used
  rawseries[rawseries[,"PopulationUntransformed"]==0, "PopulationUntransformed"] = NA
  
  # there might be missing obs: SeriesStep will tell us
  series = matrix(NA, max(rawseries[,"SeriesStep"]) + 1)
  
  # copy data at the proper time spot
  series[rawseries[, "SeriesStep"] + 1] = rawseries[, "PopulationUntransformed"]
  
  # remove NAs at begining of series
  while (is.na(series[1])) series = series[-1]
  
  # log-tranform series
  series = log(series)
  
  # vector of booleans indicating if obs is NA or not
  indices = as.integer(is.na(series))
  
  # replace NAs by 0 before feeding to admb
  series[is.na(series)] = 0
  
  list_length = c(list_length, length(series))     
  list_series = c(list_series, list(series))
  list_indices = c(list_indices, list(indices))
}

# Build matrix of observations
maxLength = max(unlist(list_length))
nbseries = length(list_length)

obs = matrix(0, nbseries, maxLength)
nas = matrix(1, nbseries, maxLength)
for (i in 1:nbseries) {
  obs[i,1:list_length[[i]]] = list_series[[i]]
  nas[i,1:list_length[[i]]] = list_indices[[i]]
}
lengths = unlist(list_length)

# OK obs, nas and lengths will do the job

# load taxonomic data
taxon = read.csv("Taxonomics.csv")
taxon[378,"TaxonomicFamily"] = "Noctuidae"

# Analysis by order
# we keep the orders with more than 10 time series
tsbyorder = aggregate(taxon$ID, by=list(taxon$TaxonomicOrder), length)

#selected = tsbyorder[order(tsbyorder$x, decreasing = T)[c(1, 3:16)], ]
selected = tsbyorder[order(tsbyorder$x, decreasing = T)[1:14], ]

order.sel = taxon$TaxonomicOrder %in% selected[,1]

# subset data
obs=obs[order.sel,]
nas=nas[order.sel,]
lengths = lengths[order.sel]

# build factors for each model
order.fact = droplevels(taxon[order.sel,"TaxonomicOrder"])
class.fact = droplevels(taxon[order.sel,"TaxonomicClass"])
family.fact = droplevels(taxon[order.sel, "TaxonomicFamily"])
full.fact = factor(rownames(taxon[order.sel,]))
commonb.fact = gl(1,sum(order.sel), labels="commonb")

save.image("TaxonomyData.Rdata")