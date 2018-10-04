# Managing Projects with Uncertain Deadlines
Tom M Logan  
www.tomlogan.co.nz

## Files

    |- main.R
        |- generate_projects.R
      
## Steps

1. main
    i. set seed
    ii. vary parameters
        * link density
        * variance of activity times
        * variance of deadline
        * baseline deadline (?)
    ii. generate projects (calls `generate_projects.R`)
    iii. compare crash strategies
    iv. record result
    v. repeat
2. generate projects
    i. generate a binary upper triangular matrix using a parameter for link density
    i. remove paths which are subpaths of other paths
    i. draw activity mean (random uniform variable (0,1)) and variance (exponential random variable (4))
3. compare strategies
    i. loop crash strategies
    ii. update activities accordingly
4. determine project properties (call `CalculatePathProperties` in `main.R`)
    i. identify paths from precedence matrix
    ii. loop through paths
        * calculate mean
        * calculate variance
        * determine covariance
5. on-time completion probability (call `CalculateCompletionProbability` in `main.R`)


## Crashing rules:
1. Crash the critical path
Calculate a crash reserve (the total amount of crashing we can use): 20% of sum of activity means
    a. Identify the critical path
		b. Identify activities on the critical path
		c. Rank the activities based on
			i. Their mean time
		d. Crash the activity
			i. Reduce mean by 50% or the remaining crashing reserve.
			ii. Variance remains constant
Repeat until crashing reserve is depleted.

2. Probability of lateness
  a. Calculate the probability that each path is completed before the deadline
  b. For each activity, calculate the product of the probabilities of its paths completing on-time
  c. For each activity, subtract this probability from 1. This gives the probability of this activity being on at least one late path.
  c. Rank by this probability, and if equal, by the mean
  d. Crash the activity
    i. Reduce the mean 50%
    ii. Variance constant
  Crash the top activities until the crashing reserve is spent.

#### Citation:

