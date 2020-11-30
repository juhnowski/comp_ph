install.packages("wesanderson")
library(ggplot2)

hamiltonian <- read.csv("~/comp_ph/6/diagonal/dia0/hamiltonian.csv")
ggplot(hamiltonian, aes(X, Y, fill= Z)) + geom_tile()

spin_mat <- read.csv("~/comp_ph/6/diagonal/dia0/spin_mat.csv")
View(spin_mat)
ggplot(spin_mat, aes(X, Y, fill= Z)) + geom_tile()

tr_mat <- read.csv("~/comp_ph/6/diagonal/dia0/tr_mat.csv")
View(tr_mat)
ggplot(tr_mat, aes(X, Y, fill= Z)) + geom_tile()

tr_vec <- read.csv("~/comp_ph/6/diagonal/dia0/tr_vec.csv")
View(tr_vec)
ggplot(tr_vec, aes(X, Y, fill= Z)) + geom_tile()

diag_mat <- read.csv("~/comp_ph/6/diagonal/dia0/diag_mat.csv")
View(diag_mat)
ggplot(diag_mat, aes(X, Y, fill= Z)) + geom_tile()

diag_vec <- read.csv("~/comp_ph/6/diagonal/dia0/diag_vec.csv")
View(diag_vec)
ggplot(diag_vec, aes(X, Y, fill= Z)) + geom_tile()

eig_s <- read.csv("~/comp_ph/6/diagonal/dia0/eig_s.csv")
View(eig_s)
ggplot(eig_s, aes(X, Y, fill= Z)) + geom_tile()

eig_m <- read.csv("~/comp_ph/6/diagonal/dia0/eig_m.csv")
View(eig_m)
ggplot(eig_m, aes(X, Y, fill= Z)) + geom_tile()