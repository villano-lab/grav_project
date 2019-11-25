These notes are from our meeting on 8 Nov, 2019. </br>
I am collecting them here so we can reference them later.</br>

Weighting in lmit:
- Any weighting method can be used
- Chi-square is calculated from sum of squares of residuals
- “Note that the calculation of chi-square and reduced chi-square assume that the returned residual function is scaled properly to the uncertainties in the data. For these statistics to be meaningful, the person writing the function to be minimized must scale them properly.”
- https://lmfit.github.io/lmfit-py/model.html residuals account for weights
- Chi-square and reduced chi-square are described here https://lmfit.github.io/lmfit-py/fitting.html#fit-results-label 
- https://stackoverflow.com/questions/43421482/python-lmfit-reduced-chi-square-too-small-after-weighted-fit states that lmfit multiplies the residual by the weight (before squaring it); this means the final equation for chi-square is consistent with the standard: http://maxwell.ucsc.edu/~drip/133/ch4.pdf -- assuming we use 1/stdev for our weights
    - I couldn’t find anything official on this part, so I will have to trust this stackexchange answer.
TL;DR if we want a standard reduced chi-square value then we should use weights equal to 1/stdev.