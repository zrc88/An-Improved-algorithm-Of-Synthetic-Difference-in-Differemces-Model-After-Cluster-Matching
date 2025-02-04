Overview

Imagine a policy implemented at the city level in China. We want to examine how the policy’s effects differ between treated and untreated cities. The traditional Difference-in-Differences (DID) approach is not ideal here because cities vary significantly and should not be treated equally. Instead, we allocate appropriate weights to the control group to construct a counterfactual for the treatment group. This new method, Synthetic DID (SDID), aims to generate a counterfactual for the treated units, yielding more precise results than traditional DID when there are significant differences at the individual or time level.

However, the original SDID method has limitations, especially when there are many treatment groups. While SDID assigns different weights to control units, it simply averages the outcomes for all treatment units—effectively giving them equal weight. Consider that the policy affects a range of cities: major metropolitan areas like Beijing and Shanghai, as well as smaller cities such as Jiamusi. Even if you have never heard of Jiamusi, you can infer that its economic conditions are very different from those of the larger cities. This discrepancy means that combining these diverse units to form a single counterfactual is not meaningful.

Proposed Method

To address this issue, I propose dividing the treatment group into smaller clusters of cities that share similar characteristics. The steps are as follows:

Clustering:
Use k-means or other clustering algorithms to group the treated cities based on relevant characteristics.

Counterfactual Construction:
For each cluster, construct an appropriate counterfactual (i.e., “anti-factor”) by allocating weights to the control units.

ATT Calculation:
Compute the Average Treatment Effect (ATT) separately for each cluster.

Aggregation:
Combine the cluster-specific ATT estimates by weighting each group according to the number of units it contains. This approach aims to produce a more precise and reasonable overall outcome than existing methods.

Computational Improvements

By using clustering algorithms, we can achieve more accurate mathematical outcomes. However, this approach requires additional computational resources. Without further optimization, the computation time may increase roughly in proportion to the number of clusters. Due to time and expertise constraints, I have employed parallel computing techniques to calculate the ATT for different clusters concurrently.

Computational Limitations

Although parallel computing helps reduce overall computation time—potentially matching or even improving upon the time cost of traditional SDID (depending on the number of available CPU cores)—the placebo tests remain the most time-consuming part of the process. At this stage, I have not been able to optimize this component further.

Methodological Limitations

Regarding the new algorithm applied after cluster matching, I believe that unit fixed effects play a much more critical role than time fixed effects in this model. However, in the current approach, I have not accounted for the heterogeneity of time weights after treatment; I have only justified the method mathematically at the individual level. I anticipate that future algorithmic advancements may better address and optimize the treatment of time effects.
