library(openxlsx) 
library(dplyr)
library(ggplot2)

#######################################################
#######################################################
### Build a function to do this
# - body.part, xyz.axis and repetition: to identify specific groups of patients
# - use.last.value: set it to TRUE to replace NA's with the last available value, 
#                   FALSE to keep the NA's.
load.patients.df <- function(filename,
                             body.part = "Knee",
                             xyz.axis = "X",
                             repetition = 1,
                             use.last.value)
{
  # Load the data frame
  patients.df <- read.csv(filename, sep = ';', )
  
  print(patients.df %>% dim())
  
  # Transpose matrix to put patients on the rows
  patients.df <- as.data.frame(t(patients.df))
  
  print(patients.df %>% dim())
  
  # First row is not meaningful
  patients.df <- patients.df[-1,]
  
  print(patients.df %>% dim())
  
  # Cut the columns not containing meaningful info
  patients.df <- patients.df %>% select(!c(V2,V3))
  
  print(patients.df %>% dim())
  
  # Take only data relative to knee
  patients.df <- patients.df %>% filter(V1 == body.part)
  
  print(patients.df %>% dim())
  
  # Take only data relative to X axis
  patients.df <- patients.df %>% filter(V4 == xyz.axis)
  
  print(patients.df %>% dim())
  
  # Keep only first repetition
  patients.df <- patients.df[grep(pattern = paste("_0", repetition, sep=""), 
                                  x = rownames(patients.df)),]
  
  print(patients.df %>% dim())
  
  # Give Names to the columns
  col.names <- c("Body.part", 
                 "Axis.angle", 
                 paste("T.", rep(1:(ncol(patients.df)-2)), sep=""))
  patients.df <- setNames(patients.df, col.names)
  
  # Remove columns not containing data (optional)
  patients.df <- patients.df %>% select(-c(Body.part, Axis.angle))
  
  print(patients.df %>% dim())
  
  # Compute length of each time series of each patient
  series_length <- rep(0, nrow(patients.df))
  for(i in 1:nrow(patients.df))
  {
    first.nan <- which.max(is.nan(as.numeric(patients.df[i,])))
    series_length[i] <- first.nan-1
  }
  
  # # Remove patients with time series shorter than 500 time steps
  # if(length(which(series_length<500)) > 0)
  # {
  #   patients.df <- patients.df[-which(series_length<500),]
  #   series_length <- series_length[-which(series_length<500)]
  # }
  
  # Replace NaN with last available datum
  if(use.last.value)
  {
    for(i in 1:nrow(patients.df))
    {
      first.nan <- which.max(is.nan(as.numeric(patients.df[i,])))
      patients.df[i, first.nan:ncol(patients.df)] <- patients.df[i, first.nan-1]
    }
  }
  
  return(list(data = patients.df, n = series_length))
}


# Function to merge datasets
merge.patients.df <- function(list.df, max.n.col, group.labels, replace.w.nan)
{
  # The merged df to return
  patients.df <- NULL

  id <- 0
  for(df in list.df)
  {
    id <- id + 1
    # Add the missing columns at the end
    curr.ncol <- ncol(df)
    
    print(nrow(df))
    
    # How many NaN to add
    to.add <- max.n.col - curr.ncol
    add.df <- NULL
      
    # Add NaN at the end
    if(replace.w.nan)
    {
      # Create columns of NaN
      add.df <- as.data.frame(matrix(NaN, nrow = nrow(df), ncol = to.add))
  
    } else # Add last available value
    {
      # At the end of the row add the last available values
      for(i in 1:nrow(df))
      {
        row.to.add <- rep(df[i,curr.ncol], to.add)
        add.df <- rbind(add.df, row.to.add)
      }
      add.df <- as.data.frame(add.df)
    }
      
    # Introduce names for the added columns
    new.names <- NULL
    if(curr.ncol < max.n.col) new.names <- paste("T.", (curr.ncol+1):max.n.col, sep="")
    else new.names <- NULL
    
    add.df <- setNames(add.df, new.names)
    
    temp.df <- NULL
    # Merge the additional data to the existing df
    if(to.add != 0)
      temp.df <- data.frame(df, add.df)
    else
      temp.df <- df
    
    # Make data numeric
    temp.df <- apply(X = temp.df, MARGIN = 2, FUN = as.numeric)
    
    # Add factor indicating the group
    temp.df <- data.frame(group.id = as.factor(rep(group.labels[id], nrow(df))), temp.df)

    # Merge the modified df into the global one
    patients.df <- as.data.frame(rbind(patients.df, temp.df))
  }
  
  return(patients.df)
}












