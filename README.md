# Managing Projects with Uncertain Deadlines
Tom M Logan  
www.tomlogan.co.nz

## Background

## Files

    |- main.R
        |- generate_projects.R
        |- compare_strategies.R
            |- determine_project_properties.R
            |- calculate_completion_probability.R
                |- calculate_true_completion.R


## Steps

1. main
    i. set seed
    ii. vary parameters
        * link density
        * variance of activity times
        * variance of deadline
        * baseline deadline (?)
    ii. generate projects (calls `generate_projects.R`)
    iii. compare crash strategies (call `compare_strategies.R`)
    iv. record result
    v. repeat
2. generate projects
    i. generate matrix paths and activities (a 1 indicates that an activity, column, is in the path, row)
    i. find all paths in the directed acyclic graph
    i. remove paths which are subpaths of other paths
    i. draw activity mean (random uniform variable (0,1)) and variance (exponential random variable (4))
2. generate projects (B)
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
    *  
    * this is either the `calc_acpa` function from the modified PERT paper or the true multivariate normal completion time
    * $P(max(G_1, .., G_n) \leq 0)$ where $G_i$ is the difference between the completion time of path $i$ and the deadline.




#### Citation:
