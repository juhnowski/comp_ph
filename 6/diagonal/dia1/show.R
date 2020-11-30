library(ggplot2)

eig_0_1 <- read.csv("~/comp_ph/6/diagonal/dia1/eig_0_1.csv")
View(eig_0_1)
ggplot(eig_0_1, aes(X, Y, fill= Z)) + geom_tile()

eig_1_8 <- read.csv("~/comp_ph/6/diagonal/dia1/eig_1_8.csv")
View(eig_1_8)
ggplot(eig_1_8, aes(X, Y, fill= Z)) + geom_tile()

eig_2_28 <- read.csv("~/comp_ph/6/diagonal/dia1/eig_2_28.csv")
View(eig_2_28)
ggplot(eig_2_28, aes(X, Y, fill= Z)) + geom_tile()

eig_3_56 <- read.csv("~/comp_ph/6/diagonal/dia1/eig_3_56.csv")
View(eig_3_56)
ggplot(eig_3_56, aes(X, Y, fill= Z)) + geom_tile()

eig_4_70 <- read.csv("~/comp_ph/6/diagonal/dia1/eig_4_70.csv")
View(eig_4_70)
ggplot(eig_4_70, aes(X, Y, fill= Z)) + geom_tile()