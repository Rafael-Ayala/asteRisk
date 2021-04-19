library(deSolve)
library(profvis)
library(asteRisk)
library(plotly)
t = 0
# Y = matrix(c(-7136143.23821097, -415951.978574370, -505103.379544365, -575.578554821520, 1079.95441989431, 7358.26105845415), nrow=6, ncol = 6)
Y = c(-7136143.23821097, -415951.978574370, -505103.379544365, -575.578554821520, 1079.95441989431, 7358.26105845415)
MJD_UTC = 52388.9135185187
solarArea = 110.5
satelliteMass = 8000
satelliteArea = 62.500000000000000
Cr = 1
Cd = 2.200000000000000

parameters = list(
    MJD_UTC = 52388.9135185187,
    solarArea = 110.5,
    satelliteMass = 8000,
    satelliteArea = 62.500000000000000,
    Cr = 1,
    Cd = 2.200000000000000
)

asteRisk:::accel(10000, Y, MJD_UTC, solarArea, satelliteMass, satelliteArea, Cr, Cd)

times <- c(0,1)

test_ode <- profvis(ode(y=Y, times=times, func=asteRisk:::odeModel, parms=parameters, method="ode45"))

times <- seq(0, by=60, length.out = 1588+1)
# test_ode2 <- profvis(ode(y=Y, times=times, func=asteRisk:::odeModel, parms=parameters, method="radau", maxsteps=1, rtol=1e-13, atol=1e-16))
test_ode2 <- profvis(ode(y=Y, times=times, func=asteRisk:::odeModel, parms=parameters, method="radau", maxsteps=1, rtol=1e-13, atol=1e-16, hini=0.01))

test_radau <- ode(y=Y, times=times, func=asteRisk:::odeModel, parms=parameters, method="radau", maxsteps=100, rtol=1e-13, atol=1e-16, hini=0.01)

#test_ode2 <- profvis(ode(y=Y, times=times, func=asteRisk:::odeModel, parms=parameters))

x <- test_radau[,2]
y <- test_radau[,3]
z <- test_radau[,4]
c <- test_radau[,1]
data <- data.frame(x,y,z,c)

fig <- plot_ly(data, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'lines+markers',
               line = list(width = 6, color = ~c, colorscale = 'Viridis'),
               marker = list(size = 3.5, color = ~c, colorscale = 'Greens', cmin = -20, cmax = 50))

fig

