Imagine there is a policy implemented in China and it acts on city level. We
want to observe how the effects are in the treated city than untreated city. DID
is not a good option to tackle this because cities are so different then should
not be considered equally. So we allocate appropriate weights to the control
groups and fit anti-factor situation of the treatment group. The new method is
Synthetic DID(SDID). SDID aims to constitute an anti-factor of the treatment
group to obtain a more precise outcome than traditional DID when groups are
relatively different in individual or time level.
However, original SDID is also limited in many situation, especially there
are many treatment groups because SDID allocates different weights to control
units, but for treatment groups, they are simply take the mean value, in other
words, the same weights. Let’s go back to the assumption. The policy affect
many different cities, including the metropolitan Beijing, Shanghai, and Jiamusi
is also included. You maybe never hear about the city but you could easily tell
it could not as rich as the two metropolitan or you have heard it before. That is
where the problem lies. You can not use it together with the two to constitute
the so-called anti-factor otherwise there is no meaning in it because the antifactor can not be defined as any city’s anti-factor.
So what I want to do is divide the treatment group into smaller groups
and each of them shares similar characteristics. I intend to use k-means(or
other algorithms) to distinguish them. After that, fit anti-factor situation for
each group and calculate the average treatment effect(ATT) dividedly. Finally,
allocate group weights based on the number of units in the group to get a more
precise and reasonable outcome than existed methods.


Improvement in Computing
By cluster algorithm we can obtain a more precious outcome in mathematics.
However, it costs more computation resources. If I do not optimize it in some
aspects, time consumption will increase by as many multiples as the number
of groups. Due to limited time and professional knowledge, I just use parallel
computing technology here when calculate different groups ATT.
Shortcoming in Computing
Just as mentioned last section, I can only simplify the calculation process
to several times or equal time cost than traditional SDID based on the cores in
CPU. But there is a remained problem that placebo computation is the most
time-consuming process in it. I can not optimize it in this stage.
Shortcoming in this method
Talking about the new algorithm in the SDID after cluster matching. Personally, I think units fixed effects are far more important than time fixed effects
in this model. However, the reality is I did not consider the time weights heterogeneity after treatment. I just justify it in mathematics and at individual level.
If there is a better algorithm is launched in the future, I think it will optimize
the time part.
