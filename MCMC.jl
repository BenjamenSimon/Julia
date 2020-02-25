using Distributions
using InvertedIndices
using Random
using Plots
using DataFrames

using Distances
using LinearAlgebra


function MH_accept(llh, llh_prime)
  alpha =  llh_prime - llh
  return( min(0, alpha) )
end


function MH_accept_MRW(llh, llh_prime, nu_prior, lambda_prior, param_cur, param_prime)

  alpha = llh_prime + (nu_prior)*log(param_prime) - lambda_prior*param_prime -
                  ( llh + (nu_prior)*log(param_cur) - lambda_prior*param_cur)

  return( min(0, alpha) )
end


function Folded_draw(param_cur, param_sigma, lower, upper)

  p_draw = rand(Normal(param_cur, param_sigma))

  while p_draw > upper || p_draw < lower

    if p_draw > upper
      p_draw = p_draw - 2*(p_draw - upper)
    end

    if p_draw < lower
       p_draw = p_draw + 2*(lower - p_draw)
    end

  end

  return(p_draw)
end


function pmatform(dist_mat, ps, ds)

  pmat = Array{Float64}(undef, size(dist_mat))

  index = findall(dist_mat .< ds)
  pmat[index] .= ps[1]

  index = findall(dist_mat .>= ds)
  pmat[index] .= ps[2]

  pmat[diagind(pmat)] .= 0

  return pmat
end


function Initialise_2b(N, rem_times = rem_times,
                         beta1_true = missing, beta2_true = missing, p_true = missing,
                         dist_true = missing, inf_times = missing,
                         d_lower = d_lower, d_upper = d_upper, dist_mat = dist_mat,
                         reparam = true)

  beta1_cur, beta2_cur, p_cur, dist_cur, beta_mat, inf_times_cur = fill(missing, 6)
  arb = true
  while arb == true
    # Initialise parameters
    if ismissing(beta1_true)
      beta1_cur = rand(Uniform(0, 0.4))
    else
      beta1_cur = beta1_true
    end
    if ismissing(beta2_true)
       beta2_cur = rand(Uniform(0, beta1_cur))
     else
       beta2_cur = beta2_true
     end
    if ismissing(p_true)
      p_cur = rand(Uniform(0, 1))
    else
      p_cur = p_true
    end
    if ismissing(dist_true)
      dist_cur = rand(Uniform(d_lower, d_upper))
    else
      dist_cur = dist_true
    end

    # Initialise the beta matrix
    if reparam == true
      beta_mat = betamatform(dist_mat, [beta1_cur, (p_cur*beta1_cur)], dist_cur)
    else
      beta_mat = betamatform(dist_mat, [beta1_cur, beta2_cur], dist_cur)
    end

    # Initialise the infection times
    if ismissing(inf_times)
      inf_periods = rand(Exponential(1/0.15), N)
      inf_times_cur = rem_times - inf_periods
    else
      inf_times_cur = inf_times
    end

    #Calculate the log-likelihood
    if simple_log_likelihood(inf_times_cur, rem_times, beta_mat) != -Inf
      break
    end
  end

  println("initialised")

  return(beta1_cur, beta2_cur, p_cur, dist_cur, inf_times_cur)
end


function MCMC_reparam(N_its, N, inf_ids, rem_times, dist_mat,
                 d_upper, sigmap, sigmad; infupdate = 1,
                 lambda_b1 = 0.001, nu_b1 = 1, lambda_g = 0.001 , nu_g = 1,
                 inc_beta1 = [missing, true], inc_p = [missing, true], inc_dist = [missing, true],
                 inc_inf_times = [missing, true], inc_gamma = [missing, true])

  ############################################
  ### Calculate upper and lower bound on d ###
  ############################################

      inf_dist_mat = deepcopy(dist_mat[inf_ids, inf_ids])
      d_lower = minimum(inf_dist_mat[inf_dist_mat .> 0])

  ##############################################
  ### Initialise the parameters and epidemic ###
  ##############################################

      InitialiseEpi = Initialise_2b(N, rem_times,
                                     inc_beta1[1], missing, inc_p[1],
                                     inc_dist[1], inc_inf_times[1],
                                     d_lower, d_upper, dist_mat,
                                     true)
      beta1_cur, beta2_cur, p_cur, dist_cur, inf_times = InitialiseEpi

      beta_mat = betamatform(dist_mat, [beta1_cur, (p_cur*beta1_cur)], dist_cur)

  ######################
  ### Results matrix ###
  ######################

      res = Array{Float64}(undef, N_its, 7)

  ##########################
  ### Functional objects ###
  ##########################

      it = 1
      n_I = length(inf_ids)
      acc_sum_I = acc_sum_p = acc_sum_d = 0
      llh = simple_log_likelihood(inf_times, rem_times, beta_mat)
      inf_inds1 = findall(inf_times .< Inf)

  ###########################
  ### ~~ THE ALGORITHM ~~ ###
  ###########################

      while it <= N_its

        ###############################
        ### Gibbs sampler for gamma ###
        ###############################

        if inc_gamma[2] == true

            #Calculate the removal integral
            g_ints = sum(rem_times[inf_ids] - inf_times[inf_ids])

            # Draw gamma
            gamma_cur = rand( Gamma( (n_I+nu_g), 1/(lambda_g + g_ints) ) )
            #rgamma(n=1, shape = (n.I+nu.g), rate = (lambda.g + g.ints))

        end


        ###############################
        ### Gibbs sampler for beta1 ###
        ###############################

        if inc_beta1[2] == true

            # Form the p matrix
            p_mat = pmatform(dist_mat, [1, p_cur], dist_cur)

            # Calculate the infection integral
            b_ints = simple_integral_part(inf_times, rem_times, p_mat)


            # Draw beta1
            beta1_cur = rand( Gamma( (n_I-1+nu_b1), 1/(lambda_b1 + b_ints) ) )
            #<- rgamma(n=1, shape = (n.I-1+nu.b1), rate = (lambda.b1 + b.ints))

            # Calculate functional objects
            beta_mat = betamatform(dist_mat, [beta1_cur, (p_cur*beta1_cur)], dist_cur)
            llh = simple_log_likelihood(inf_times, rem_times, beta_mat)

        end



        ###################################
        ### MH step for infection times ###
        ###################################

        if inc_inf_times[2] == true

            # Which infection time is being replaced
            Ireplace = sample(inf_ids, infupdate, replace = false)

            # Draw new infection time
            Qdraw = rand(Exponential(1/gamma_cur), infupdate)
            inf_times_prime = deepcopy(inf_times)
            inf_times_prime[Ireplace] = deepcopy(rem_times[Ireplace] - Qdraw)

            # Calculate functional objects
            llh_prime = simple_log_likelihood(inf_times_prime, rem_times, beta_mat)

            # MH acceptance probability
            alpha_I = MH_accept(llh, llh_prime)

            # Do we accept the new time(s) or not?
            accept_test = log(rand())
            if accept_test < alpha_I  #If yes:

              inf_times = deepcopy(inf_times_prime)
              llh = deepcopy(llh_prime)

              acc_sum_I = acc_sum_I+1
            end
        end


        #####################
        ### MH step for p ###
        #####################

        if inc_p[2] == true

            # Draw new p value
            p_draw = Folded_draw(p_cur, sigmap, 0, 1)

            # Calculate functional objects
            beta_mat_prime = betamatform(dist_mat, [beta1_cur, (p_draw * beta1_cur)], dist_cur)
            llh_prime = simple_log_likelihood(inf_times, rem_times, beta_mat_prime)

            # MH acceptance probability
            alpha_p = MH_accept(llh, llh_prime)

            # Do we accept the new p or not?
            accept_test = log(rand())
            if accept_test < alpha_p  #If yes:

              p_cur = deepcopy(p_draw)
              beta_mat = deepcopy(beta_mat_prime)
              llh = deepcopy(llh_prime)

              acc_sum_p = acc_sum_p+1
          end
        end

        #####################
        ### MH step for d ###
        #####################

        if inc_dist[2] == true

            # Draw new d value
            dist_draw = Folded_draw(dist_cur, sigmad, d_lower, d_upper)

            # Calculate functional objects
            beta_mat_prime = betamatform(dist_mat, [beta1_cur, (p_cur*beta1_cur)], dist_draw)
            llh_prime = simple_log_likelihood(inf_times, rem_times, beta_mat_prime)

            # MH acceptance probability
            alpha_d = MH_accept(llh, llh_prime)

            # Do we accept the new d or not?
            accept_test = log(rand())
            if accept_test < alpha_d  #If yes:

              dist_cur = deepcopy(dist_draw)
              beta_mat = deepcopy(beta_mat_prime)
              llh = deepcopy(llh_prime)

              acc_sum_d = acc_sum_d+1
            end
        end

        ##########################
        ### Record the results ###
        ##########################

            # Record parameters
            res[it,:] = [it, beta1_cur, (p_cur * beta1_cur), gamma_cur, dist_cur, p_cur , llh]

            # Update count
            it = it + 1
      end

  res = DataFrame(res)
  names!(res, [:sample, :beta1, :beta2, :gamma, :d, :p, :llh])

  return(res, (acc_sum_I/N_its), (acc_sum_p/N_its), (acc_sum_d/N_its))
end









#####################
####### TEST ########
#####################


#seedfinder(100, [0.004, 0.002], 2, 25, 2, 20)

Random.seed!(3)
exN = 100
exdistmat = unifdistmat(exN, 20, 20)[2]
exbetamat = betamatform(exdistmat, [0.004 0.002], 2)
exresults = GSEsim(exN, exbetamat, 0.15)

exinf_ids = findall(exresults[:,2] .< Inf)
exinf_times = exresults[:,2]
exrem_times = exresults[:,3]


res1, accI, accp, accd = MCMC_reparam(100000, 100, exinf_ids, exrem_times, exdistmat,
                 15, 0.1, 1, infupdate = 1,
                 lambda_b1 = 0.001, nu_b1 = 1, lambda_g = 0.001 , nu_g = 1,
                 inc_beta1 = [missing, true], inc_p = [missing, true], inc_dist = [missing, true],
                 inc_inf_times = [exinf_times, true], inc_gamma = [missing, true])

plot(res1[:,2], title = "beta 1")
plot(res1[:,3], title = "beta 2")
plot(res1[:,4], title = "gamma")
plot(res1[:,5], title = "d")
plot(res1[:,6], title = "p")
plot(res1[:,7], title = "llh")

plot(res1[:,2], res1[:,5], title = "d")

histogram(res1[:,2], title = "beta 1")
histogram(res1[:,3], title = "beta 2")
histogram(res1[:,4], title = "gamma_cur")
histogram(res1[:,5], title = "d")
histogram(res1[:,6], title = "p")
histogram(res1[:,7], title = "llh")

describe(res1[5001:100000, :])

@time begin
  MCMC_reparam(10000, 100, exinf_ids, exrem_times, exdistmat,
                   15, 0.1, 1, infupdate = 3,
                   lambda_b1 = 0.001, nu_b1 = 1, lambda_g = 0.001 , nu_g = 1,
                   inc_beta1 = [missing, true], inc_p = [missing, true], inc_dist = [missing, true],
                   inc_inf_times = [exinf_times, false], inc_gamma = [missing, true])
end
