# library(dslabs)
# data(heights)
# heights
# nrow(heights)
# heights[1032,2]
# which.min(heights$height)
# median(heights$height)
# x <- heights %>% group_by(sex) %>% summarise(n())
# sum(heights$height > 78 & heights$sex == "Female")
library(tidyverse)
library(caret)
library(dslabs)
data(heights)
y <- heights$sex
x <- heights$height
set.seed(2007)
test_index <- createDataPartition(y,times = 1, list = FALSE)
test_set <- heights[test_index,]
train_set <- heights[-test_index]
