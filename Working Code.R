
# To do: CALIBRATION using Lux SMC

# Done: 1. Initialized LOB
      # 2. Simulated LOB for 250 time steps
      # 3. Calculated the mid-price at the end of each time step
      # 4. midprice will be used for callibration


# Steps taken to achieve this:
      # A) SIMULATED PRICE of LOs:
          # 1. q_taker (Mean reverting random walk - depends on delta_s)
          # 2. lambda(t) (depends on q_taker (MC simulations and current), C_lam and lambda0)
          # 3. nu(t) (exponetial lambda(t) R.V.)
          # 4. The Limit Order Book (Buy and Ask side with pb and pa = best bid and ask, respectively)

      # B) SIMULATED QUANTITY of LOs and MOs 
        # 1. Simulate a number of LIMIT orders in each time period (depending on N_a and alpha)
        # 2. Then, Initialize the order book for 10 MC steps with only LOs hitting the LOB
        # 3. Execute MOs

      # C) Simulated Entire Order Book for 250 time steps


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

lambda0 = 0.7 # initial value of lambda
C_lam = 0.15

lambda_generate = function(lambda0, C_lam, mean_q, q_taker0, delta_s, Num_orders){
  lambda = lambda0
  q_taker = q_taker_evolve(TimeHorizon = Num_orders, mean_q, q_taker0, delta_s)

  for (k in 1:((Num_orders+1))){
    mean_MC_q = Monte_Carlo_q_taker(num_steps=10^5, TimeHorizon=1, mean_q, q_taker0, delta_s)
    lambda[k+1] = lambda0*(1 + abs(q_taker[k+1]-0.5)/sqrt(mean_MC_q) * C_lam)
  }
  return(list(lambda, q_taker))
}

lambda_q_taker = lambda_generate(lambda0, C_lam, mean_q, q_taker0, delta_s, Num_orders=10)
plot(lambda_q_taker[[1]], type="l")

# now we want to simulate nu: exponentially distributed random number

nu_for_price = function(TimeHorizon, lambda){
  nu = c()
  for (k in 1:(TimeHorizon+1)){
    uniform = runif(n=1, 0,1)
    nu[k] = -lambda[k]*log(uniform)
  }
  return(nu)
}

nu = nu_for_price(TimeHorizon, lambda_q_taker[[1]])
plot(nu, type="l")

###################################################
# Create a limit buy and limit sell order book
# create a limit buy order book

library(dplyr)
create_LOB_0 = function(){ # requires dplyr
  p0 = 100; pb = p0; pa = p0

  Limit_Buy_init = as.data.frame(c(pb))
  colnames(Limit_Buy_init) = "price"

  Limit_Buy_Book = Limit_Buy_init%>%
    group_by(price) %>%
    summarise(count = n())

  Limit_Ask_init = as.data.frame(c(pa))
  colnames(Limit_Ask_init) = "price"

  Limit_Ask_Book = Limit_Ask_init%>%
    group_by(price) %>%
    summarise(count = n())

  myLO_Books = list(Limit_Buy_Book, Limit_Buy_Book)
  
  return(myLO_Books)
}
LOB_0 = create_LOB_0()
#######################################

# function to determine the price of the limit order 
# NB! given the current LO book with corresponding best ask (bid) price (for limit buy (sell) orders)
# and given current nu: defined nu_k

q_provider = 0.5
nu_k = nu_for_price(1, lambda_q_taker[[1]])

update_LOB_aftereach_LO = function(q_provider, myLO_Books, nu_k){
  Limit_Buy_Book = myLO_Books[[1]]
  Limit_Ask_Book = myLO_Books[[2]]
  uniform = runif(1, 0, 1)
  if (uniform < q_provider){ # we have a buy order
    LO_price = round(Limit_Ask_Book[1,1][[1]] - 1 - nu_k, 2) # round to 2 cents
    LOs = rep(Limit_Buy_Book$price, Limit_Buy_Book$count)
    Limit_Buy_Book = as.data.frame(c(LOs, LO_price))
    colnames(Limit_Buy_Book) = "price"
    Limit_Buy_Book = Limit_Buy_Book%>%
      group_by(price) %>%
      summarise(count = n()) 
  } else{ # we have a limit order
    LO_price = round(Limit_Buy_Book[nrow(Limit_Buy_Book),1][[1]] + 1 + nu_k, 2) # round to cents
    LOs = rep(Limit_Ask_Book$price, Limit_Ask_Book$count)
    Limit_Ask_Book = as.data.frame(c(LOs, LO_price))
    colnames(Limit_Ask_Book) = "price"
    Limit_Ask_Book = Limit_Ask_Book%>%
      group_by(price) %>%
      summarise(count = n()) 
  }
  return(list(Limit_Buy_Book, Limit_Ask_Book))
}

myLO_Books = update_LOB_aftereach_LO(q_provider, myLO_Books = LOB_0, nu_k= 0.5)

#################################################################################

# LIMIT ORDER FREQUENCYs

Na = 250 # Number of liquidity providers
alpha = 0.15 # frequency of MOs or LOs

Num_orders = floor(Na*alpha)

#################################################################################

# Need to create/initialize the LOB using 10 MC steps where only LOs can be placed

init_LOB = function(Na, alpha, LOB_0, q_provider, lambda0, C_lam, mean_q, q_taker0, delta_s){
  # initialize: Fill the order book before the commencement of actual trading.
  # generate and insert initial limit orders for 10 
  # Monte Carlo steps
  # The number of orders inserted depends on the number of traders, N_A, and frequency of limit orders, alpha, stored in simRun.
  # After the completion of initialization, the modified PreisModel
  # object is returned.
  
  # Generate Initial Orders for 10 MC steps
  num_LOs = 10*floor(Na*alpha) # 10 x number of LOs in each MC step
  
  # Generate lambda to obtain nu for each LO
  lambda_q_taker = lambda_generate(lambda0, C_lam, mean_q, q_taker0,delta_s, Num_orders = num_LOs) # SLOW
  lambda = lambda_q_taker[[1]]
  nu_k = nu_for_price(TimeHorizon = num_LOs, lambda)
  
  # For each LO, update LOB using nu and best ask and best bid
  myLO_Books = LOB_0
  for (i in 1:num_LOs){
    myLO_Books = update_LOB_aftereach_LO(q_provider, myLO_Books, nu_k = nu_k[i])
  }
  return(myLO_Books)
}

LOB_0 = create_LOB_0()
myLO_Books = init_LOB(Na, alpha, LOB_0, q_provider, lambda0, C_lam, mean_q, q_taker0, delta_s)
myLOB_init = myLO_Books
########################################################################

# Next we need a function that executes a single MO
# for each MC step there are [mu*Na] MOs
# probability MO = buy order = q_taker (which is not fixed - changes each MC step
# However, qtaker is fixed for all MOs in each time step --> qtaker(t)

q_taker_k = q_taker[10]
# Function to execute a single MO
execute_single_MO = function(q_taker_k, myLO_Books){
  u = runif(1, min=0, max=1)
  if (u < q_taker_k){ # market buy order hits best ask
    myLO_Books[[2]][1,2] = myLO_Books[[2]][1,2] - 1 # subtract 1 from count of best ask
    myLO_Books[[2]] = myLO_Books[[2]][myLO_Books[[2]]$count !=0, ] # get rid of rows with 0 count
  }else{ # market sell order hits best bid
    myLO_Books[[1]][nrow(myLO_Books[[1]]),2] = myLO_Books[[1]][nrow(myLO_Books[[1]]),2] - 1 # subtract 1 from count of best ask
    myLO_Books[[1]] = myLO_Books[[1]][myLO_Books[[1]]$count !=0, ] # get rid of rows with 0 count
  }
  return(myLO_Books)
}

mu = 0.0067
# function to execute all MOs in a single time step --> [Na*mu] MOs
execute_MOs_1timestep = function(myLO_Books, q_taker_k, Na, mu){
  num_MOs = floor(Na*mu)
  for (k in 1:num_MOs){
    myLO_Books = execute_single_MO(q_taker_k, myLO_Books)
  }
  return(myLO_Books)
}

myLO_Books= execute_MOs_1timestep(myLO_Books, q_taker_k, Na, mu)

####################################################################

# Now I have a function which executes MOs each time step and
# a function which executes single LOs

# So, we need a function which excutes [Na*alpha] LOs FIRST
# then the MOs 
# in each timestep
# and then returns the new LOB

LO_MO_1timestep = function(Na, alpha,mu, myLO_Books, q_provider,
                           lambda_t, nu_t, q_taker_t){
  # 1. LOs hit LOB
  num_LOs = floor(Na*alpha) # number of LOs
  for (k in 1:num_LOs){ # update LOB for each LO in the time step
    myLO_Books = update_LOB_aftereach_LO(q_provider, myLO_Books, nu_t)
  }
  
  # 2. MOs hit LOB
  myLO_Books = execute_MOs_1timestep(myLO_Books, q_taker_t, Na, mu)
  
  return(myLO_Books)
}



######################################################################

# Now that we have an algorithm which executes all LOs and then MOs in a SINGLE TIME STEP
# we need a function that executes all LOs and MO
# and returns the mid-price which is the time series that will be used to calibrate the model

Sim_OB = function(lambda0, C_lam, mean_q, q_taker0, delta_s, TimeHorizon, Na, alpha, mu, myLOB_init, q_provider){
  lambda_q_taker = lambda_generate(lambda0, C_lam, mean_q, q_taker0, delta_s, Num_orders = TimeHorizon)
  lambda = lambda_q_taker[[1]]
  q_taker = lambda_q_taker[[2]]
  nu = nu_for_price(TimeHorizon = TimeHorizon, lambda)
  
  myLO_Books = myLOB_init
  midprice = c()
  
  for (t in 1:TimeHorizon){
    lambda_t = lambda[t]; nu_t = nu[t]; q_taker_t = q_taker[t]
    myLO_Books = LO_MO_1timestep(Na, alpha, mu, myLO_Books, q_provider,
                                            lambda_t, nu_t, q_taker_t)
    midprice[t] = 0.5 * (myLO_Books[[1]]$price[nrow(myLO_Books[[1]])]+myLO_Books[[2]]$price[1]) # calculate the mid-price each time step
  }
  
  return(list(myLO_Books, midprice))
}


simLOB = Sim_OB(lambda0, C_lam, mean_q, q_taker0, delta_s, TimeHorizon=250, Na, alpha, mu, myLOB_init, q_provider)
sim.midprice = simLOB[[2]]
plot(sim.midprice, type = "l")
simLOB[[1]]
####################################################################

#