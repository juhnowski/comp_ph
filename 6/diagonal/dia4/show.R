library(ggplot2)

data <- read.csv("~/comp_ph/6/diagonal/dia4/eig_0_-1_-1_2.csv")
View(data)
ggplot(data, aes(X, Y, fill= Z)) + geom_tile()

data <- read.csv("~/comp_ph/6/diagonal/dia4/eig_0_1_-1_1.csv")
View(data)
ggplot(data, aes(X, Y, fill= Z)) + geom_tile()

data <- read.csv("~/comp_ph/6/diagonal/dia4/eig_0_1_1_7.csv")
View(data)
ggplot(data, aes(X, Y, fill= Z)) + geom_tile()

data <- read.csv("~/comp_ph/6/diagonal/dia4/eig_1_1_-1_5.csv")
View(data)
ggplot(data, aes(X, Y, fill= Z)) + geom_tile()

data <- read.csv("~/comp_ph/6/diagonal/dia4/eig_1_1_1_3.csv")
View(data)
ggplot(data, aes(X, Y, fill= Z)) + geom_tile()

data <- read.csv("~/comp_ph/6/diagonal/dia4/eig_2_1_-1_4.csv")
View(data)
ggplot(data, aes(X, Y, fill= Z)) + geom_tile()

data <- read.csv("~/comp_ph/6/diagonal/dia4/eig_2_1_1_5.csv")
View(data)
ggplot(data, aes(X, Y, fill= Z)) + geom_tile()

data <- read.csv("~/comp_ph/6/diagonal/dia4/eig_3_1_-1_5.csv")
View(data)
ggplot(data, aes(X, Y, fill= Z)) + geom_tile()

data <- read.csv("~/comp_ph/6/diagonal/dia4/eig_3_1_1_3.csv")
View(data)
ggplot(data, aes(X, Y, fill= Z)) + geom_tile()

data <- read.csv("~/comp_ph/6/diagonal/dia4/eig_4_-1_-1_4.csv")
View(data)
ggplot(data, aes(X, Y, fill= Z)) + geom_tile()

data <- read.csv("~/comp_ph/6/diagonal/dia4/eig_4_-1_1_1.csv")
View(data)
ggplot(data, aes(X, Y, fill= Z)) + geom_tile()

data <- read.csv("~/comp_ph/6/diagonal/dia4/eig_4_1_1_5.csv")
View(data)
ggplot(data, aes(X, Y, fill= Z)) + geom_tile()
