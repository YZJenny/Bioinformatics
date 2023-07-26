library('gsDesign')

control_median_DFS=9.5
intervention_meidan_DFS=13

accrual=12
followup=12*2
droupout=0.1


f <- nSurv(
  sided = 1,
  alpha = 0.2/2, #onde-side
  beta = 1-0.8,
  
  ratio=1,
  lambdaC = log(control_median_DFS)/intervention_meidan_DFS,
  hr=control_median_DFS/intervention_meidan_DFS,
  ## above for events
  
  T=accrual+followup,
  minfup = followup,
  eta = -log(1-droupout)/12 #per month
  
)
f



library('gsDesign')

control_median_DFS <- 9.5
intervention_median_DFS <- 13
accrual <- 12 # month
followup <- 12 * 2
droupout <- 0.1

# Calculate the desired hazard ratio based on the condition that intervention group's median DFS is at least 9.5 months longer
desired_hr <- (control_median_DFS + 9.5) / intervention_median_DFS

f <- nSurv(
  sided = 1,
  alpha = 0.2/2, # one-side
  beta = 1 - 0.8,
  
  ratio = 1,
  lambdaC = log(control_median_DFS) / intervention_median_DFS,
  hr = desired_hr,  # Use the desired hazard ratio for effectiveness testing
  ## above for events
  
  T = accrual + followup,
  minfup = followup,
  eta = -log(1-droupout)/12
)
f

