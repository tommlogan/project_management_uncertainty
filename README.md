# Managing Projects with Uncertain Deadlines
Tom M Logan  
www.tomlogan.co.nz

#### Please cite this code as:  
R.F. Bordley et al., Managing projects with uncertain deadlines, European Journal of Operational Research
(2018), https://doi.org/10.1016/j.ejor.2018.09.036

## Files

    |- main.R
        |- generate_projects.R
      
## Steps

1. main  
    1. set seed  
    2. vary parameters  
        * link density  
        * variance of activity times  
        * variance of deadline  
        * baseline deadline (?)  
    3. generate projects (calls `generate_projects.R`)  
    4. compare crash strategies  
    5. record result  
    6. repeat
2. generate projects
    1. generate a binary upper triangular matrix using a parameter for link density
    2. remove paths which are subpaths of other paths
    3. draw activity mean (random uniform variable (0,1)) and variance (exponential random variable (4))
3. compare strategies
    1. loop crash strategies
    2. update activities accordingly
4. determine project properties (call `CalculatePathProperties` in `main.R`)
    1. identify paths from precedence matrix
    2. loop through paths
        * calculate mean
        * calculate variance
        * determine covariance
5. on-time completion probability (call `CalculateCompletionProbability` in `main.R`)


## Crashing rules:
* Crash the critical path  
Calculate a crash reserve (the total amount of crashing we can use): 20% of sum of activity means  
    1. Identify the critical path  
    2. Identify activities on the critical path  
    3. Rank the activities based on  
        * Their mean time  
    5. Crash the activity    
        1. Reduce mean by 50% or the remaining crashing reserve.  
        2. Variance remains constant  
Repeat until crashing reserve is depleted.

* Probability of lateness:  
    1. Calculate the probability that each path is completed before the deadline  
    2. For each activity, calculate the product of the probabilities of its paths completing on-time  
    3. For each activity, subtract this probability from 1. This gives the probability of this activity being on at least one late path.  
    4. Rank by this probability, and if equal, by the mean  
    5. Crash the activity    
        1. Reduce the mean 50%  
        2. Variance constant  
  Crash the top activities until the crashing reserve is spent.


