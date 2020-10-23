#  http://www.sthda.com/english/articles/36-classification-methods-essentials/146-discriminant-analysis-essentials-in-r/
#
#  Ejemplo de discriminante que encontre
#
#  species sería el modulo al que pertenece cada trophospecies y tendriamos dos variables (equivalentes a Sepal.Length etc...) el nivel trofico que es continua y Habitat que es discreta y no podria usarse con el discriminante, pero parece ser que lo usan igual.
#
# https://stats.stackexchange.com/questions/158772/can-we-use-categorical-independent-variable-in-discriminant-analysis
#
# el otro tema es que no encontré como se hace la permutacion para testear, supestamente con el paquete vegan 


install.packages('caret')
library(tidyverse)
library(caret)
theme_set(theme_classic())
data("iris")
# Split the data into training (80%) and test set (20%)
set.seed(123)
training.samples <- iris$Species %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data <- iris[training.samples, ]
test.data <- iris[-training.samples, ]


# Estimate preprocessing parameters
preproc.param <- train.data %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train.data)
test.transformed <- preproc.param %>% predict(test.data)

library(MASS)
# Fit the model ----------------------------> species sería el modulo al que pertenece cada y lo que se 
model <- lda(Species~., data = train.transformed)
model
plot(model)


predictions <- model %>% predict(test.transformed)
names(predictions)

lda.data <- cbind(train.transformed, predict(model)$x)
ggplot(lda.data, aes(LD1, LD2)) +
  geom_point(aes(color = Species))


