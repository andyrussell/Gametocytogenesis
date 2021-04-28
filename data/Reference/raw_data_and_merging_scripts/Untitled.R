## Script for rendering reference files

## 







## DOZI-binders dataset

## from here: https://datascienceplus.com/extracting-tables-from-pdfs-in-r-using-the-tabulizer-package/ 
## 
library(tabulizer)
library(dplyr)

# Location of WARN notice pdf file
location <- 'https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-014-0493-0/MediaObjects/13059_2014_493_MOESM1_ESM.pdf'

# Extract the table
out <- extract_tables(location)

final <- data.frame()
for (i in 1:length(out)) {            # 2. sequence
  if  (dim(out[[i]])[2] == 17) {
    df <- data.frame(out[[i]])[ , c(1, 8, 9, 10)]
    names(df) <- c('id', 'CITH', 'DOZI', 'CD')
    final <- rbind(final, )
  } else if (dim(out[[i]])[2] == 14) {
    df <- data.frame(out[[i]])[ , c(1, 8, 9, 10)]
    names(df) <- c('id', 'CITH', 'DOZI', 'CD')
    final <- rbind(final, )
  } else if  (dim(out[[i]])[2] == 16) {
    df <- data.frame(out[[i]])[ , c(1, 7, 8, 9)]
    names(df) <- c('id', 'CITH', 'DOZI', 'CD')
    final <- rbind(final, df)
  } else if  (dim(out[[i]])[2] == 13) {
    df <- data.frame(out[[i]])[ , c(1, 7, 8, 9)]
    names(df) <- c('id', 'CITH', 'DOZI', 'CD')
    final <- rbind(final, df)
  } else if  (dim(out[[i]])[2] == 15) {
    df <- data.frame(out[[i]])[ , c(1, 6, 7, 8)]
    names(df) <- c('id', 'CITH', 'DOZI', 'CD')
    final <- rbind(final, df)
  } else if  (dim(out[[i]])[2] == 12) {
    df <-  data.frame(out[[i]])[ , c(1, 6, 7, 8)]
    names(df) <- c('id', 'CITH', 'DOZI', 'CD')
    final <- rbind(final, df)
  }
}

final <- do.call(rbind, out[-length(out)])
library(data.table)
final <- rbindlist(out, fill = TRUE)


# table headers get extracted as rows with bad formatting. Dump them.
final <- as.data.frame(final[3:nrow(final), ])

# Column names
#headers <- c('Notice.Date', 'Effective.Date', 'Received.Date', 'Company', 'City', 'No.of.Employees', 'Layoff/Closure')

# Apply custom column names
names(final) <- headers

## inspect
head(final)