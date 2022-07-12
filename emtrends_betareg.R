
require("stats")
require("glmmTMB")
require("betareg")
require("rstan")
require("rstanarm")
require("emmeans")
require("ggplot2")
require("ggrepel")

x = rnorm(100, 0, 1)
e = rnorm(100, 0, 1)
B = c(0, 3.14)
p = c(0.1, 0.1)

mu = binomial(link = logit)$linkinv(B[1] + B[2]*x)
phi = binomial(link = log)$linkinv(p[1] + p[2]*e)
y = rbeta(100, mu*phi, (1-mu)*phi)

dd = data.frame(x = x, y = (y*(length(y)-1) + 0.5)/length(y) )

m.glm = glmmTMB(y~x, data = dd, family = beta_family)
m.betareg = betareg(y~x, data = dd)
m.rstan = rstanarm::stan_betareg(y ~ x, data = dd)

slopes.1 = c(
  glm=fixef(m.glm)[[1]][2],
  betareg=coef(m.betareg)[2],
  rstan=coef(m.rstan)[2]
)

slopes.2 = c(
  glm=as.data.frame(emtrends(m.glm, ~1, var = "x"))[1,2],
  betareg=as.data.frame(emtrends(m.betareg, ~1, var = "x"))[1,2],
  rstan=as.data.frame(emtrends(m.rstan, ~1, var = "x"))[1,2]
)

gdf = data.frame(
  model.output = slopes.1,
  emtrend = slopes.2,
  model = names(slopes.1)
)

ggplot(gdf, aes(x = model.output, y = emtrend, label = model)) +
  geom_point() +
  geom_label_repel() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  xlim(c(0,2)) +
  ylim(c(0,2)) +
  theme_bw()
