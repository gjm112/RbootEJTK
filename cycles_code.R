x <- seq(0,48,1)
#y <- cos(2*pi*x/12) + rnorm(length(x),0,10)
a <- 12
s <- 12
p <- 24

#y <- (-1/asymmetry)*((x - phase_shift) %% phase + 1)*(((x - phase_shift) %% phase) < asymmetry) + ((1/(phase - asymmetry))*((x %% phase) - phase_shift) - (asymmetry/(phase - asymmetry)))*(((x %% phase) - phase_shift) > (asymmetry))
#y <- (-1/a * x + 1)*(x < a) + (1/(p-a) * x - (a/(p-a)))*(x > a)
y <- (-1/a * (x-s) %% p + 1)*((x-s) %% p < a) + (1/(p-a) * (x-s) %% p - (a/(p-a)))*((x-s) %% p > a)

#y <- cos(2*pi*(x - phase_shift)/(2*asymmetry)) * ((x - phase_shift) <= asymmetry) + 
  #cos(2*pi*(x - phase_shift - phase)/(2*(phase - asymmetry))) * ((x - phase_shift) > asymmetry) + 
  #rnorm(length(x), 0, 0)

x <- seq(0,48,1) 
plot(x, y, type = "l")
abline(v = a)
abline(h = 0)