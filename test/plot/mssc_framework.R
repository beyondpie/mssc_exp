## plot mssc framework
lwd <- 8

x <- seq(-4, 4, length = 100)
y <- dnorm(x)
plot(x, y, type = "l", axes = FALSE, lwd = lwd, xlab = "", ylab = "")

xa <- seq(-4, 4, length = 100)
ya <- dnorm(xa, mean = -0.05)
plot(xa, ya, type = "l", axes = FALSE, lwd = lwd, xlab = "", ylab = "", col = "red")

xb <- seq(-4, 4, length = 100)
yb <- dnorm(xb, mean = 0.05)
lines(xb, yb, type = "l", lwd = lwd, col = "blue")


## case 1
x1 <- seq(-4, 4, length = 100)
y1 <- dnorm(x = x1, mean = -1.2, sd = 1, log = FALSE)
plot(x, y, type = "l", axes = FALSE, lwd = lwd, xlab = "", ylab = "", col = "#F7931E" )

x2 <- seq(-4, 4, length = 100)
y2 <- dnorm(x = x2, mean = -0.7, sd = 1, log = FALSE)
lines(x2, y2, type = "l", lwd = lwd, col = "#F15A24")

x3 <- seq(-4, 4, length = 100)
y3 <- dnorm(x = x3, mean = -1, sd = 1, log = FALSE)
lines(x3, y3, type = "l", lwd = lwd, col = "#C1272D")

## case 2
x4 <- seq(-4,4, length = 100)
y4 <- dnorm(x = x4, mean = 1.2, sd = 1, log = FALSE)
lines(x4, y4, type= "l", lwd = lwd, col = "#00A99D")

x5 <- seq(-4, 4, length = 100)
y5 <- dnorm(x = x5, mean = 1, sd = 1, log = FALSE)
lines(x5, y5, type = "l", lwd = lwd, col = "#29ABE2")

x6 <- seq(-4, 4, length = 100)
y6 <- dnorm(x = x6, mean = 0.7, sd = 1, log = FALSE)
lines(x6, y6, type = "l", lwd = lwd, col = "#2E3192")


## batch effect
x1 <- seq(-4, 4, length = 100)
y1 <- dnorm(x = x1, mean = -1, sd = 1, log = FALSE)
plot(x, y, type = "l", axes = FALSE, lwd = lwd, xlab = "", ylab = "", col = "#F7931E" )

x2 <- seq(-4, 4, length = 100)
y2 <- dnorm(x = x2, mean = -0.7, sd = 1, log = FALSE)
lines(x2, y2, type = "l", lwd = lwd, col = "#F15A24")

x3 <- seq(-4, 4, length = 100)
y3 <- dnorm(x = x3, mean = -0.4, sd = 1, log = FALSE)
lines(x3, y3, type = "l", lwd = lwd, col = "#C1272D")

x4 <- seq(-4,4, length = 100)
y4 <- dnorm(x = x4, mean = 1, sd = 1, log = FALSE)
lines(x4, y4, type= "l", lwd = lwd, col = "#00A99D")

x5 <- seq(-4, 4, length = 100)
y5 <- dnorm(x = x5, mean = 0.5, sd = 1, log = FALSE)
lines(x5, y5, type = "l", lwd = lwd, col = "#29ABE2")

x6 <- seq(-4, 4, length = 100)
y6 <- dnorm(x = x6, mean = 0.7, sd = 1, log = FALSE)
lines(x6, y6, type = "l", lwd = lwd, col = "#2E3192")

