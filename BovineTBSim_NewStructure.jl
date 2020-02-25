using Distributions
using InvertedIndices
using Random
using DataFrames

################################
### Functions to create farm ###
################################

function create_farm_res(parish, farm, Time, cS, cI, bS, bI)

  farm_res = fill(0., Time, 10)

  #[:parish, :farm, :t, :cS, :cE, :cI, :bS, :bE, :bI, :move_res]

  farm_res[1, :] .= [parish, farm, 1., cS, 0., cI, bS, 0., bI, false]

  return(farm_res)
end

#farm_res = create_farm_res(2, 3, 10, 50, 2,30, 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function create_farm_track(Time)

   #[:t,
   #:sus_on, :exp_on, :inf_on, :sus_off, :exp_off, :inf_off,
   #:env, :c_exp_prob, :b_exp_prob, :c_inf_prob, :b_inf_prob,
   #:c_new_exp, :c_new_inf, :b_new_exp, :b_new_inf,
   #:E_detected, :I_detected,
   #:c_birth, :c_S_death, :c_E_death, :c_I_death,
   #:b_birth, :b_S_death, :b_E_death, :b_I_death,
   #:test_date]

  farm_track = fill(0., Time, 27)
  farm_track[1,:] .= [1.; fill(0., 25); sample(1:(365*2), 1)]

  return(farm_track)
end

#farm_track = create_farm_track(10)

##################################
### Functions to create parish ###
##################################


function create_parish_res(parish, Time, farms)

  # Inputs
  #[:parish, :farm, :Time, :cS, :cI, :bS, :bI]
  inputs = fill(0., farms, 7)

  inputs[:, 1] .= deepcopy(parish)
  inputs[:, 3] .= Time
  inputs[:, 2] = 1:farms
  inputs[:, 4] = sample(50:200, farms, replace = true)
  inputs[:, 5] = rand(Binomial(1, 0.1), farms)
  inputs[:, 6] = sample(10:50, farms, replace = true)
  inputs[:, 7] = rand(Binomial(1, 0.1), farms)

  # Parish

  parish_res = Array{Array{Float64, 2}}(undef, 1, farms)

  for i = 1:farms
    parish_res[1, i] = create_farm_res(inputs[i, 1],
      inputs[i, 2],
      convert(Int64, inputs[i, 3]),
      inputs[i, 4],
      inputs[i, 5],
      inputs[i, 6],
      inputs[i, 7])
  end

  return(parish_res)
end

#parish_res = create_parish_res(2., 10, 5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function create_parish_track(Time, farms)

  parish_track = Array{Array{Float64, 2}}(undef, 1, farms)

  # Parish
  for i = 1:farms
    parish_track[1, i] = create_farm_track(convert(Int64, Time))
  end

  return(parish_track)
end

#parish_track = create_parish_track(10, 5)

###################################
### Functions to create country ###
###################################

function create_country_res(Time, parishes)

  # Inputs
  #[:parish, :Time, :farms]
  inputs = fill(0., parishes, 3)

  inputs[:, 1] = 1:parishes
  inputs[:, 2] .= Time
  inputs[:, 3] = sample(5:15, parishes, replace = true)


  # Country
  country_res = Array{Array{Array{Float64,2}, 2}}(undef, 1, parishes)

  for i = 1:parishes
    country_res[1, i] = create_parish_res(inputs[i, 1],
                                          convert(Int64, inputs[i, 2]),
                                          convert(Int64, inputs[i, 3]))
  end
  return(country_res)
end

#country_res = create_country_res(10, 3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function create_country_track(Time, parishes, farms)

  # inputs
    inputs = fill(0., parishes, 2)

    inputs[:,1] .= Time
    inputs[:,2] = farms

    # Country
    country_track = Array{Array{Array{Float64,2}, 2}}(undef, 1, parishes)

    for i = 1:parishes
      country_track[1, i] = create_parish_track(inputs[i, 1], convert(Int64, inputs[i, 2]))
    end
    return(country_track)
end

#n_parishes = size(country_res, 2)
#n_farms = fill(0, n_parishes)
#for i in 1:n_parishes
#  global n_farms[i]  = size(country_res[1, i], 2)
#end

#country_track = create_country_track(10, 3, n_farms)




#######################################
### Functions to generate movements ###
#######################################

function extract_states(t, n_parishes, n_farms)
  #[:parish, :farm, :t, :cS, :cE, :cI, :bS, :bE, :bI, :move_res]
  states = fill(0., sum(n_farms), 10)

  index = 1
  for i = 1:n_parishes
    for j = 1:n_farms[i]
      for k = 1:10
        begin
          states[index, k] = country_res[1, i][1, j][t, k]
          states
        end
      end
      index = index + 1
    end
  end
  return (states)
end

#current_states = extract_states(1, n_parishes, n_farms)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function rng_mvhyper(n, k)
  # n is a vector of size of pop for each group
  # k is the number of trials (total number moving)
  N = sum(n) # total pop size
  m = length(n) # number of groups
  n_otr = N - n[1] # number not in first group

  x = fill(0, m) # results

  x[1] = rand(Hypergeometric(n[1], n_otr, k))

  for i in 2:(m-1)
    n_otr = n_otr - n[i]
    k = k - x[i-1]
    x[i] = rand(Hypergeometric(n[i], n_otr, k))
  end

  x[m] = k - x[m-1]
  return(x)
end

#rng_mvhyper([150, 10, 20], 10)

function generate_moves(t, current_states)
  #[:parish, :farm, :t, :cS, :cE, :cI,
  #  :cN, :cMove_off, :cS_move_off, :cE_move_off, :cI_move_off,
  #  :cS_move_on, :cE_move_on, :cI_move_on, :move_res]
  movements = fill(0., sum(n_farms), 15)

  # :parish, :farm, :t, :cS, :cE, :cI
  movements[:,1:6] .= current_states[:, 1:6]
  # :cN
  movements[:, 7] .= sum.(eachrow(current_states[:,4:6]))

  index = 1
  for i in movements[:, 7] # :cN
    # :cMove_off
    movements[index, 8] = rand(Binomial(convert(Int64, i), (0.05/31)), 1)[1]
    index = index+1
  end

  index2 = 1
  for j in movements[:, 8] # :cMove_off
    # :cS_move_off, :cE_move_off, :cI_move_off
    movements[index2, 9:11] .= rng_mvhyper(movements[index2, 4:6], j)
    index2 = index2 + 1
  end

  # movement restrictions
  movements[:, 15] .= current_states[:, 10]

  move_free = findall(movements[:, 15] .== 0)
  movements[Not(move_free), 8:11] .= 0

  # :cS_move_on, :cE_move_on, :cI_move_on
  movements[:, 12:14] .= 0

  return(movements)
end

#movements = generate_moves(1, current_states)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function assign_moves(movements)
  move_free = findall(movements[:, 15] .== 0)
  moved = findall(movements[move_free, 8] .> 0)
  moved = move_free[moved]

  for i in moved
    parish = convert(Int64, movements[i, 1])
    new_parish = sample((1:n_parishes)[Not(parish)])

    farm = sample(findall(movements[move_free,1] .== new_parish))
    farm = move_free[farm]

    movements[farm, 12] = movements[farm, 12] + movements[i, 9]
    movements[farm, 13] = movements[farm, 13] + movements[i, 10]
    movements[farm, 14] = movements[farm, 14] + movements[i, 11]
  end

  return(movements)
end

#movements = assign_moves(movements)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

function moves_t(t, current_states)
  gen_moves = generate_moves(t, current_states)
  moves = assign_moves(gen_moves)

  return(moves)
end

#movements = moves_t(1)

##############################################################
### Functions for creating and tracking parish environment ###
##############################################################

function create_p_env(Time)
  n_parishes = size(country_res, 2)

  p_env = fill(0., Time, n_parishes)

  return(p_env)
end

#p_env = create_p_env(10)

function calculate_p_env(t, parish, ep, current_states)
  env = p_env[t, parish]

  parish_mems = findall(current_states[:,1] .== parish)

  sumI = sum(sum.(eachcol(current_states[parish_mems, [6,9]])))
  sumN = sum(sum.(eachcol(current_states[parish_mems, 4:9])))

  env = env*exp(-ep) + 1/ep * sumI/sumN * (1 - exp(-ep))

  return(env)
end

#calculate_p_env(1, 3, 0.05)


#########################################
### Function to generate timestep t+1 ###
#########################################

function update_states(t, parish, farm, p_env, movement)

  #####################
  ### Find the farm ###
  #####################

  farms = 1:sum(n_farms)

  farm_i_res = deepcopy(country_res[1, parish][1, farm])

  farm_i_track = deepcopy(country_track[1, parish][1, farm])

  move_row = findall(x -> movement[x,1] == parish && movement[x,2]== farm, farms)


  ##############################
  ### Initialise the records ###
  ##############################

  new_state = deepcopy(farm_i_res[t, :])
  new_state[Not(1:2)] .= 0

  new_track = deepcopy(farm_i_track[t, :])
  new_track[Not(1)] .= 0

  ##############################
  ### Extract current states ###
  ##############################

  c_S = farm_i_res[t, 4]
  c_E = farm_i_res[t, 5]
  c_I = farm_i_res[t, 6]
  c_N = c_S + c_E + c_I

  b_S = farm_i_res[t, 7]
  b_E = farm_i_res[t, 8]
  b_I = farm_i_res[t, 9]
  b_N = b_S + b_E + b_I

  test_date = farm_i_track[t, 27]
  move_res = farm_i_res[t, 10]

  #########################################
  ### Calculate environmental reservoir ###
  #########################################

  env = farm_i_track[t, 8]
  env = env*exp(-ep) + 1/ep * (c_I+b_I)/(c_N+b_N) * (1 - exp(-ep))

  parish_env = p_env[t, parish]

  ########################
  ### Update movements ###
  ########################

  sus_move_off = (movement[move_row, 9])[1]
  exp_move_off = (movement[move_row, 10])[1]
  inf_move_off = (movement[move_row, 11])[1]

  sus_move_on = (movement[move_row, 12])[1]
  exp_move_on = (movement[move_row, 13])[1]
  inf_move_on = (movement[move_row, 14])[1]

  c_S_new = c_S + sus_move_on - sus_move_off
  c_E_new = c_E + exp_move_on - exp_move_off
  c_I_new = c_I + inf_move_on - inf_move_off
  c_N_new = c_S_new + c_E_new + c_I_new

  ######################################################
  ### Calculate exposure and infection probabilities ###
  ######################################################

  if c_N_new > 0
    c_exp_prob = 1 - exp( - c_I_new/c_N_new * c_beta - f * env - pf * parish_env)
  else
    c_exp_prob = 0
  end

  b_exp_prob = 1 - exp( - b_I/b_N * b_beta - f * env - pf * parish_env)

  c_inf_prob = 1 - exp( - gamma)
  b_inf_prob = 1 - exp( - gamma)

  #####################################################
  ### Draw new exposures and infectious individuals ###
  #####################################################

  c_new_exp = rand(Binomial(convert(Int64, c_S_new), c_exp_prob))
  b_new_exp = rand(Binomial(convert(Int64, b_S), b_exp_prob))

  c_new_inf = rand(Binomial(convert(Int64, c_E_new), c_inf_prob))
  b_new_inf = rand(Binomial(convert(Int64, b_E), b_inf_prob))

  c_S_new = c_S_new - c_new_exp
  c_E_new = c_E_new + c_new_exp - c_new_inf
  c_I_new = c_I_new + c_new_inf
  c_N_new = c_S_new + c_E_new + c_I_new

  b_S_new = b_S - b_new_exp
  b_E_new = b_E + b_new_exp - b_new_inf
  b_I_new = b_I + b_new_inf
  b_N_new = b_S_new + b_E_new + b_I_new

  #######################
  ### Draw detections ###
  #######################
  global I_detected
  global E_detected
  if (I_detected + E_detected) > 0
    E_detected = 0
    I_detected = 0
  end

  if t == test_date
    E_detected = rand(Binomial(convert(Int64, c_E_new), edet))
    I_detected = rand(Binomial(convert(Int64, c_I_new), det))

    if (I_detected + E_detected) > 0
      test_date = test_date + 30
      move_res = 1 #true
    else
      test_date = test_date + 365*2
      move_res = 0 #false
    end
  end

  c_E_new = c_E_new - E_detected
  c_I_new = c_I_new - I_detected
  c_N_new = c_S_new + c_E_new + c_I_new

  ##############################
  ### Draw births and deaths ###
  ##############################

  c_birth = rand(Poisson(c_N_new*0.25/365))

  c_S_death = rand(Binomial(convert(Int64, c_S_new), (0.25/365)))
  c_E_death = rand(Binomial(convert(Int64, c_E_new), (0.25/365)))
  c_I_death = rand(Binomial(convert(Int64, c_I_new), (0.25/365)))

  b_birth = rand(Poisson(b_N_new*0.25/365))

  b_S_death = rand(Binomial(convert(Int64, b_S_new), (0.25/365)))
  b_E_death = rand(Binomial(convert(Int64, b_E_new), (0.25/365)))
  b_I_death = rand(Binomial(convert(Int64, b_I_new), (0.25/365)))

  c_S_new = c_S_new - c_S_death + c_birth
  c_E_new = c_E_new - c_E_death
  c_I_new = c_I_new - c_I_death

  b_S_new = b_S_new - b_S_death + b_birth
  b_E_new = b_E_new - b_E_death
  b_I_new = b_I_new - b_I_death

  #########################
  ### Record new states ###
  #########################

  new_state[3] = t+1
  new_track[1] = t+1

  new_state[4] = c_S_new
  new_state[5] = c_E_new
  new_state[6] = c_I_new

  new_state[7] = b_S_new
  new_state[8] = b_E_new
  new_state[9] = b_I_new
  new_state[10] = move_res

  ##########################################
  ### Record useful tracking information ###
  ##########################################

  new_track[Not(1)] .= [sus_move_on, exp_move_on, inf_move_on, sus_move_off, exp_move_off, inf_move_off,
                       env, c_exp_prob, b_exp_prob, c_inf_prob, b_inf_prob,
                       c_new_exp, c_new_inf, b_new_exp, b_new_inf,
                       E_detected, I_detected,
                       c_birth, c_S_death, c_E_death, c_I_death,
                       b_birth, b_S_death, b_E_death, b_I_death, test_date]

  return([new_state, new_track])
end

#update = update_states(1, 1, 3, p_env, movements)


###################
### The Process ###
###################

## The Parameters ##

c_beta = 0.00002
b_beta = 0.00004
gamma = 0.015

f = 0.0004
pf = 0.0001
ep = 0.05

det = 0.75
edet = 0.2*det

## Initialise the country ##

Time = 365*4+1
n_parishes = 10

country_res = create_country_res(Time, n_parishes)

n_farms = fill(0, n_parishes)
for i in 1:n_parishes
  global n_farms[i]  = size(country_res[1, i], 2)
end

country_track = create_country_track(Time, n_parishes, n_farms)

## Functional objects ##

#farms_vec <- farms_per_parish %>% map(~1:.x) %>% unlist()
#Farm_IDs <- data.frame(parish = rep(1:parishes, times = farms_per_parish), farm = farms_vec, nrow = 1:sum(farms_per_parish))

env = 0
p_env = create_p_env(Time)

E_detected = 0
I_detected = 0

t = 1



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

while t < Time

  global current_states_top = extract_states(t, n_parishes, n_farms)

  global movements_top = moves_t(t, current_states_top)

  for i in 1:n_parishes
    for j in 1:n_farms[i]
      update = update_states(t, i, j, p_env, movements_top)
        for k in 1:10
          country_res[1,i][1,j][(t+1),k] = update[1][k]
          country_track[1,i][1,j][(t+1),k] = update[2][k]
        end
    end
  end

  for l in 1:n_parishes
    p_env[(t+1), l] = calculate_p_env(t, l, ep, current_states_top)
  end

  global t = t+1

end

end_states = extract_states(t, n_parishes, n_farms)
