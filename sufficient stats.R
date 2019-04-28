# Algorithm to find sufficient statistics

# 1. s() = I()
# 2. Compute s(y_obs) = y_obs
y_obs = runif(200, min=95, max = 105)
# 3. Sample N thetas from pi(theta)
# Prior dbn:
#   delta, tilda, mu U(0,0.1)
#   alpha U(0.1, 0.5)
#   lambda_0 U(0,200)
#   C_lam  U(0,20)

N = 50
num = N
prior_sim = function(num){
  deltas = runif(num, 0, 0.1)
  tildas = runif(num, 0, 0.1)
  mus = runif(num, 0, 0.1)
  alphas = runif(num, 0.1, 0.5)
  lambda_0s = runif(num, 0, 2)
  C_lams = runif(num, 0, 1)
  
  prior_mat = cbind(deltas, tildas, mus, alphas, lambda_0s, C_lams)
  
  return(prior_mat)
}

prior_mat = prior_sim(N)
head(prior_mat)

# 4. Using priors, simulate observations
N = 50
yobs = matrix(NA, nrow = N, ncol = TimeHorizon)
for (m in 1:N){
  simLOB = Sim_OB(lambda0, C_lam, mean_q, q_taker0, delta_s, TimeHorizon=250, Na, alpha, mu, myLOB_init, q_provider)
  yobs[m,] = simLOB[[2]]
}

plot(yobs[5,], type = "l")
# 4. Regress theta = b0 + b1y1 + b2y2 + ... + bTYT
reg_df = cbind(prior_mat[1:N,1], yobs)
fit = lm(prior_mat[,1]~ , data =)