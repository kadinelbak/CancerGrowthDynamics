# Models module - Contains all ODE model definitions
module Models

export logistic_growth!, logistic_growth_with_death!, gompertz_growth!, 
       gompertz_growth_with_death!, exponential_growth_with_delay!, 
       logistic_growth_with_delay!, exponential_growth!, 
       exponential_growth_with_death_and_delay!

# 1) plain logistic: p = (r, K)
function logistic_growth!(du,u,p,t)
  r,K = p
  K = max(K, eps(eltype(u)))
  du[1] = r*u[1]*(1 - u[1]/K)
end

# 2) logistic + death: p = (r, K, death_rate)
function logistic_growth_with_death!(du,u,p,t)
  r,K,death_rate = p
  K = max(K, eps(eltype(u)))
  du[1] = r*u[1]*(1 - u[1]/K) - death_rate*u[1]
end

# 3) Gompertz: p = (a, b, K) - safe implementation
function gompertz_growth!(du,u,p,t)
  a,b,K = p
  if u[1] <= 0 || u[1] >= K
    du[1] = 0.0
  else
    du[1] = a*u[1]*log(K/u[1])
  end
end

# 4) Gompertz + death: p = (a, b, K, death_rate) - safe implementation  
function gompertz_growth_with_death!(du,u,p,t)
  a,b,K,death_rate = p
  if u[1] <= 0 || u[1] >= K
    du[1] = -death_rate*u[1]
  else
    du[1] = a*u[1]*log(K/u[1]) - death_rate*u[1]
  end
end

# 5) exp with delay: p = (r, K, t_lag) - corrected to include carrying capacity
function exponential_growth_with_delay!(du,u,p,t)
  r,K,tlag = p
  K = max(K, eps(eltype(u)))
  du[1] = (t>=tlag ? r : 0.0)*u[1]*(1 - u[1]/K)
end

# 6) logistic with delay: p = (r, K, t_lag)
function logistic_growth_with_delay!(du,u,p,t)
  r,K,tlag = p
  K = max(K, eps(eltype(u)))
  du[1] = (t>=tlag ? r : 0.0)*u[1]*(1-u[1]/K)
end

# 7) pure exponential: p = (r,)
function exponential_growth!(du,u,p,t)
  r = p[1]; du[1] = r*u[1]
end

# 8) exponential with death and delay: p = (r, K, death_rate, t_lag)
function exponential_growth_with_death_and_delay!(du,u,p,t)
  r,K,death_rate,tlag = p
  K = max(K, eps(eltype(u)))
  du[1] = (t>=tlag ? r : 0.0)*u[1]*(1 - u[1]/K) - death_rate*u[1]
end

end # module Models
