# To do:  1. Simulate a number of LIMIT orders in each time period (depending on N_a and alpha)
        # 2. Then, Initialize the order book for 10 MC steps with only LOs hitting the LOB
        # 3. Execute MOs

# Done: Thus far, I have simulated the PRICE of limit buy and ask orders which depend on 
        # 1. q_taker (Mean reverting random walk - depends on delta_s)
        # 2. lambda(t) (depends on q_taker (MC simulations and current), C_lam and lambda0)
        # 3. nu(t) (exponetial lambda(t) R.V.)
        # 4. The Limit Order Book (Buy and Ask side with pb and pa = best bid and ask, respectively)


###################################################################
# PRICE
    # lambda(t) = current mkt depth
            # lambda0 = init placement depth param
            # C_lam = integer param
            # <[q_taker - 0.5]^2> = iterate q_taker for 10^5 MC steps
            # and take mean of [q_taker - 0.5]^2


    # nu = [-lambda(t)ln(u)], where U~U(0,1)
    # p_b = best bid
    # p_a = best ask


    # p_LB = price of limit buy order
    #      = p_a - 1 - nu
    # p_LS = price of limit sell order
    #      = p_b + 1 + nu

# NUMBER of LOs:
  # N_A = num. liquidity providers = Number of Liquidity Suppliers
  # alpha = frequency of LOs
  # for every MC step [alpha*N_A] LOs of SIZE 1
  # buy order wp q_provider = 0.5


#################################################################

# Simulating the mean reverting q_taker process

delta_s = 0.0328       # incremental change in probability of market buy order
TimeHorizon = 250    # number of MC steps
q_taker0 = 0.5      # initial market order buy probability is set to half
mean_q = 0.5        # mean with which the mean reverting random walk reverts to 

# Here the market buy order probability evolves over time.
# It is a mean reverting random walk with stepsize delta_s, mean of 0.5
# and mean reversion probability of 0.5 + |q - 0.5|
  
q_taker_evolve = function(TimeHorizon, mean_q, q_taker0, delta_s){
  q_taker_sim = q_taker0  # setting up the market order buy probability vector
  for (i in 1:TimeHorizon){
    mean_rev_prob = mean_q + abs(q_taker_sim[i] - mean_q)
    u01 = runif(n=1,0,1)
    if (q_taker_sim[i] < mean_q){                  # below the mean
      if (u01 < mean_rev_prob){
        q_taker_sim[i+1] = q_taker_sim[i]+delta_s
      }
      else{
        q_taker_sim[i+1] = q_taker_sim[i]-delta_s
      }
    } 
    else{                                   # above the mean
      if (u01 < mean_rev_prob){
        q_taker_sim[i+1] = q_taker_sim[i]-delta_s
      }
      else{
        q_taker_sim[i+1] = q_taker_sim[i]+delta_s
      }
    }
  }
  return(q_taker_sim)
}

q_taker = q_taker_evolve(TimeHorizon, mean_q, q_taker0, delta_s)

# Now we want to estimate an evolving lambda:

# For each lambda(t) we require 10^5 MC simulations of q_taker
# then take the mean of (qtaker - 0.5)^2
Monte_Carlo_q_taker = function(num_steps = 10^5, TimeHorizon, mean_q, q_taker0, delta_s){
  MC_q = c()
  for (j in 1:num_steps){
    q_taker_store = q_taker_evolve(TimeHorizon, mean_q, q_taker0, delta_s)
    MC_q[j] = (q_taker_store[TimeHorizon+1]-0.5)^2
  }
  return(mean(MC_q))
}

lambda0 = 100 # initial value of lambda
C_lam = 10
lambda = lambda0

for (k in 1:((TimeHorizon+1)){
  mean_MC_q = Monte_Carlo_q_taker(num_steps=10^5, TimeHorizon=1, mean_q, q_taker0, delta_s)
  lambda[k+1] = lambda0*(1 + abs(q_taker[k+1]-0.5)/sqrt(mean_MC_q) * C_lam)
}

plot(lambda, type="l")

# now we want to simulate nu: exponentially distributed random number

nu_for_price = function(TimeHorizon, lambda){
  nu = c()
  for (k in 1:(TimeHorizon+1)){
    uniform = runif(n=1, 0,1)
    nu[k] = -lambda[k]*log(uniform)
  }
  return(nu)
}

nu = nu_for_price(TimeHorizon, lambda)
plot(nu, type="l")

###################################################
# Create a limit buy and limit sell order book
# create a limit buy order book
p0 = 100;pb = p0; pa = p0

Limit_Buy_init = as.data.frame(c(pb))
colnames(Limit_Buy_init) = "price"

library(dplyr)
Limit_Buy_Book = Limit_Buy_init%>%
  group_by(price) %>%
  summarise(count = n())

Limit_Ask_init = as.data.frame(c(pa))
colnames(Limit_Ask_init) = "price"

library(dplyr)
Limit_Ask_Book = Limit_Ask_init%>%
  group_by(price) %>%
  summarise(count = n())
#######################################

# function to determine the price of the limit order 
# NB! given the current LO book with corresponding best ask (bid) price (for limit buy (sell) orders)
# and given current nu: defined nu_k

q_provider = 0.5

update_LOB_aftereach_LO = function(q_provider, myLO_Books, nu_k){
  Limit_Buy_Book = myLO_Books[[1]]
  Limit_Ask_Book = myLO_Books[[2]]
  uniform = runif(1, 0,1)
  if (uniform < q_provider){ # we have a buy order
    LO_price = Limit_Ask_Book[1,1][[1]] - 1 - nu_k
    LOs = rep(Limit_Buy_Book$price, Limit_Buy_Book$count)
    Limit_Buy_Book = as.data.frame(c(LOs, LO_price))
    colnames(Limit_Buy_Book) = "price"
    Limit_Buy_Book = Limit_Buy_Book%>%
      group_by(price) %>%
      summarise(count = n()) 
  } else{
    LO_price = Limit_Buy_Book[nrow(Limit_Buy_Book),1][[1]] - 1 - nu_k
    LOs = rep(Limit_Ask_Book$price, Limit_Ask_Book$count)
    Limit_Ask_Book = as.data.frame(c(LOs, LO_price))
    colnames(Limit_Ask_Book) = "price"
    Limit_Ask_Book = Limit_Ask_Book%>%
      group_by(price) %>%
      summarise(count = n()) 
  }
  return(list(Limit_Buy_Book, Limit_Ask_Book))
}

myLO_Books = update_LOB_aftereach_LO(q_provider, myLO_Books, nu_k= 0.5)

#################################################################################

# LIMIT ORDER FREQUENCYs


# Need to create the LOB using 10 MC steps where only LOs can be placed




########################################################################



